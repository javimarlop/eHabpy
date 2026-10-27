"""Export of Earth Engine images to Google Drive (batch tasks) and import of the exported files onto a grid.

A batch task has no limit of size or time per request, unlike the direct download of grid.download: it suits big
grids (e.g. an ecoregion at 100 m). The files reach the Drive of the Earth Engine account, in FOLDER, as NAME.tif or,
when the image is bigger than FILE_DIM pixels in a direction, as pieces NAME-<row offset>-<column offset>.tif. They
are read from the local copy of Google Drive for desktop (local_drives) or from a folder where they were downloaded,
and mosaicked onto the grid.
"""

import glob
import math
import os
import re

import ee
from osgeo import gdal

from . import areas
from .grid import NODATA

gdal.UseExceptions()

FILE_DIM = 16384  # side of the pieces of the exported files (pixels, a multiple of 256): 1 GB of float32
DONE, FAILED = 'SUCCEEDED', ('FAILED', 'CANCELLED', 'UNKNOWN')


def export(image, grid, name, folder, project, nodata=NODATA):
	"""Start a task that exports the first band of IMAGE on GRID to FOLDER of Google Drive as NAME (float32, NoData
	= NODATA). Returns the name of its operation (projects/<project>/operations/<id>)."""
	img = ee.Image(image).select([0]).toFloat().unmask(nodata, False)
	task = ee.batch.Export.image.toDrive(
		image=img, description=re.sub(r'[^A-Za-z0-9_,:;-]', '_', name)[:100], folder=folder, fileNamePrefix=name,
		region=grid.region(), crs=grid.crs, crsTransform=[grid.res, 0, grid.xmin, 0, -grid.res, grid.ymax],
		maxPixels=int(1e13), fileDimensions=[FILE_DIM, FILE_DIM], fileFormat='GeoTIFF')
	task.start()
	return task.operation_name or f'projects/{project}/operations/{task.id}'


def status(operation):
	"""(state, progress 0-1, error message) of a task: state PENDING, RUNNING, SUCCEEDED, FAILED, CANCELLED or
	UNKNOWN (not found)."""
	try:
		op = ee.data.getOperation(operation)
	except ee.EEException as e:
		return 'UNKNOWN', 0.0, str(e)
	meta = op.get('metadata', {})
	return meta.get('state', 'UNKNOWN'), float(meta.get('progress') or 0), (op.get('error') or {}).get('message')


def local_drives():
	"""Folders 'My Drive' of Google Drive for desktop on this computer (macOS and Windows)."""
	names = ('My Drive', 'Mi unidad')
	home = os.path.expanduser('~')
	patterns = [os.path.join(home, 'Library', 'CloudStorage', 'GoogleDrive-*', n) for n in names]
	patterns += [f'{d}:\\{n}' for d in 'DEFGHIJKLMNOPQRSTUVWXYZ' for n in names] if os.name == 'nt' else []
	return sorted({p for pat in patterns for p in glob.glob(pat) if os.path.isdir(p)})


def pieces(grid):
	"""Number of files of an export of GRID."""
	return math.ceil(grid.width / FILE_DIM) * math.ceil(grid.height / FILE_DIM)


def folders(folder):
	"""The local FOLDER and its namesakes 'FOLDER (1)', 'FOLDER (2)', ...: Earth Engine creates a new Drive folder
	for each task that starts before the folder exists, and Google Drive for desktop shows them with these names."""
	parent, base = os.path.split(os.path.normpath(folder))
	copies = [d for d in glob.glob(os.path.join(parent, glob.escape(base) + ' (*)'))
	          if re.fullmatch(re.escape(base) + r' \(\d+\)', os.path.basename(d))]
	return [folder] + sorted(copies)


def find(folder, name):
	"""Exported files of NAME (NAME.tif or its pieces NAME-*-*.tif) in the local FOLDER or its namesakes."""
	files = []
	for d in folders(folder):
		files += glob.glob(os.path.join(d, glob.escape(name) + '.tif'))
		files += [f for f in glob.glob(os.path.join(d, glob.escape(name) + '-*.tif'))
		          if re.fullmatch(re.escape(name) + r'-\d+-\d+\.tif', os.path.basename(f))]
	return sorted(files)


def mosaic(files, grid, path, nodata=NODATA, metadata=None):
	"""Write the exported FILES onto GRID as a float32 GeoTIFF at PATH (written as PATH.part.tif, then renamed)."""
	part = path + '.part.tif'
	try:
		ds = gdal.Warp(part, files, format='GTiff', dstSRS=areas.srs_of(grid.crs).ExportToWkt(),
		               outputBounds=grid.extent, xRes=grid.res, yRes=grid.res, resampleAlg='near',
		               srcNodata=nodata, dstNodata=nodata, outputType=gdal.GDT_Float32, multithread=True,
		               creationOptions=['COMPRESS=DEFLATE', 'PREDICTOR=3', 'TILED=YES', 'BIGTIFF=IF_SAFER'])
		if metadata:
			ds.SetMetadata({k: str(v) for k, v in metadata.items()})
		try:
			ds.GetRasterBand(1).ComputeStatistics(False)
		except RuntimeError:  # no valid pixels
			pass
		ds = None
		os.replace(part, path)
	except BaseException:
		ds = None
		if os.path.exists(part):
			os.remove(part)
		raise
