"""Grids aligned to the pixel size, and download of Earth Engine images on them as GeoTIFF, in tiles.

The edges of a grid are multiples of its pixel size, so grids of the same CRS and pixel size are aligned. Each tile
is computed by Earth Engine with ee.data.computePixels (the same computation as the app, without exporting to
Drive); a tile that fails because it is too big or too slow is split into four, down to MIN_TILE pixels.
"""

import concurrent.futures as cf
import math
import os
import time

import ee
import numpy as np
from osgeo import gdal

from . import areas

gdal.UseExceptions()

NODATA = -9999
MIN_TILE = 64
SPLIT_ERRORS = ('timed out', 'memory limit', 'too many concurrent', 'payload size', 'too large', 'must be less than')


class Grid:
	"""Grid of square pixels of size RES in CRS, with its upper left corner at (XMIN, YMAX)."""

	def __init__(self, crs, res, xmin, ymax, width, height):
		self.crs, self.res = crs, float(res)
		self.xmin, self.ymax = float(xmin), float(ymax)
		self.width, self.height = int(width), int(height)

	@classmethod
	def around(cls, crs, res, xmin, ymin, xmax, ymax, margin=0.0):
		"""Smallest grid with edges at multiples of RES that covers the box grown by MARGIN on every side."""
		res = float(res)
		x0 = math.floor((xmin - margin) / res) * res
		y0 = math.floor((ymin - margin) / res) * res
		x1 = math.ceil((xmax + margin) / res) * res
		y1 = math.ceil((ymax + margin) / res) * res
		return cls(crs, res, x0, y1, round((x1 - x0) / res), round((y1 - y0) / res))

	@property
	def pixels(self):
		return self.width * self.height

	@property
	def geotransform(self):
		return (self.xmin, self.res, 0.0, self.ymax, 0.0, -self.res)

	@property
	def extent(self):
		return self.xmin, self.ymax - self.height * self.res, self.xmin + self.width * self.res, self.ymax

	def describe(self):
		return {'crs': self.crs, 'resolution': self.res, 'width': self.width, 'height': self.height,
		        'extent': list(self.extent)}

	def region(self):
		"""Earth Engine rectangle of the grid (in EPSG:4326, for filterBounds)."""
		xmin, ymin, xmax, ymax = self.extent
		return ee.Geometry.Rectangle([xmin, ymin, xmax, ymax], self.crs, False)

	def ee_grid(self, col, row, width, height):
		"""Pixel grid of a tile for computePixels."""
		grid = {
			'dimensions': {'width': width, 'height': height},
			'affineTransform': {'scaleX': self.res, 'shearX': 0, 'translateX': self.xmin + col * self.res,
			                    'shearY': 0, 'scaleY': -self.res, 'translateY': self.ymax - row * self.res},
		}
		if self.crs.upper().startswith(('EPSG:', 'ESRI:', 'SR-ORG:')):
			grid['crsCode'] = self.crs
		else:
			grid['crsWkt'] = areas.srs_of(self.crs).ExportToWkt()
		return grid

	def tiles(self, size):
		return [(c, r, min(size, self.width - c), min(size, self.height - r))
		        for r in range(0, self.height, size) for c in range(0, self.width, size)]


def _split(tile):
	c, r, w, h = tile
	w1, h1 = max(w // 2, 1), max(h // 2, 1)
	parts = [(c, r, w1, h1), (c + w1, r, w - w1, h1), (c, r + h1, w1, h - h1), (c + w1, r + h1, w - w1, h - h1)]
	return [p for p in parts if p[2] > 0 and p[3] > 0]


def _compute(image, grid, tile, retries=4):
	"""Values of a tile (2D float32 array), split into four tiles if it fails because it is too big or too slow."""
	c, r, w, h = tile
	for attempt in range(retries):
		try:
			arr = ee.data.computePixels({'expression': image, 'fileFormat': 'NUMPY_NDARRAY',
			                             'grid': grid.ee_grid(c, r, w, h)})
			return arr[arr.dtype.names[0]].astype(np.float32)
		except ee.EEException as e:
			msg = str(e).lower()
			if any(s in msg for s in SPLIT_ERRORS) and max(w, h) > MIN_TILE:
				out = np.empty((h, w), np.float32)
				for p in _split(tile):
					out[p[1] - r:p[1] - r + p[3], p[0] - c:p[0] - c + p[2]] = _compute(image, grid, p, retries)
				return out
			if attempt == retries - 1 or not any(s in msg for s in ('429', 'quota', 'unavailable', '503', '500',
			                                                          'deadline', 'internal error')):
				raise
		except (OSError, ConnectionError):
			if attempt == retries - 1:
				raise
		time.sleep(2 ** attempt * 5)


def download(image, grid, path, tile=512, workers=8, nodata=NODATA, progress=None, metadata=None):
	"""Download the first band of IMAGE on GRID as a float32 GeoTIFF at PATH (NoData = NODATA), in tiles of TILE
	pixels computed by WORKERS threads. PROGRESS(done, total) is called after each tile. The file is written as
	PATH.part.tif and renamed when complete."""
	img = ee.Image(image).select([0]).toFloat().unmask(nodata, False)
	part = path + '.part.tif'
	drv = gdal.GetDriverByName('GTiff')
	ds = drv.Create(part, grid.width, grid.height, 1, gdal.GDT_Float32,
	                ['COMPRESS=DEFLATE', 'PREDICTOR=3', 'TILED=YES', 'BIGTIFF=IF_SAFER'])
	ds.SetGeoTransform(grid.geotransform)
	ds.SetProjection(areas.srs_of(grid.crs).ExportToWkt())
	band = ds.GetRasterBand(1)
	band.SetNoDataValue(nodata)
	if metadata:
		ds.SetMetadata({k: str(v) for k, v in metadata.items()})
	tiles = grid.tiles(tile)
	pool = cf.ThreadPoolExecutor(max_workers=workers)
	try:
		futures = {pool.submit(_compute, img, grid, t): t for t in tiles}
		for i, fut in enumerate(cf.as_completed(futures), 1):
			c, r, _, _ = futures[fut]
			band.WriteArray(fut.result(), c, r)
			if progress:
				progress(i, len(tiles))
		pool.shutdown()
		try:
			band.ComputeStatistics(False)
		except RuntimeError:  # no valid pixels
			pass
		ds.FlushCache()
		ds = band = None
		os.replace(part, path)
	except BaseException:
		pool.shutdown(wait=False, cancel_futures=True)
		ds = band = None
		if os.path.exists(part):
			os.remove(part)
		raise
