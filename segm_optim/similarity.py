#### eHabitat+ HFT module: similarity of the HFTs in their ecoregion and landscape metrics of the similar areas
#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2026 (from ehab_optim.py and pas/subpas_loop_segm_optim.py)
#### License: CC BY-SA 3.0

"""
Similarity of every HFT of a scenario with the rest of its ecoregion, landscape metrics of the similar areas and
the indices of the study area (thdi.R). Run from the segm_optim folder, after run_hft.py:

    python3 similarity.py                         # the options of config.py (ECOSIM_*)
    python3 similarity.py --cut p50 --res 500     # another cut of the tree, with 500 m pixels
    python3 similarity.py --plan                  # only show the grid and the pixels of every HFT (JSON)
    python3 similarity.py --area gee_snevada/ecoregion   # compare them with another area (a folder of inVars)

For every HFT, the pixels of the similarity grid (the pixel of the variables of inVars/<set>/ecoregion, or a
multiple of it with ECOSIM_RESOLUTION) that the HFT covers (at least ECOSIM_MIN_COVER of their area, from the
polygons of the HFTs of the cut) give its mean and covariance of the variables. Every pixel of the ecoregion gets
the chi-square probability of its Mahalanobis distance to the HFT (0 if it is 0.001 or less): the similarity. The
pixels outside the HFT with a similarity at least as high as the mean (or median, ECOSIM_STAT) of the HFT pixels
are the similar areas, and their patches (ECOSIM_NEIGHBOURS) give the landscape metrics of eHabitat+ (HSR,
number of patches, largest patch, aggregation...). thdi.R then computes the indices of each study area (THDI) and
similarity_plots.R the plots.

The area compared is the folder ecoregion of the variables of the HFTs (the RESOLVE ecoregions that the study area
touches, as downloaded by get_vars_gee.py) or, with ECOSIM_AREA, any folder with the same variables: a comparison
area, inside or outside the study area (get_vars_gee.py --area). The statistics of the HFTs come from its subfolder
study_area when it has one (the study area at the pixel of the area, which get_vars_gee.py downloads when the area
does not contain the study area), else from the area itself.

The results go to results/<STUDY_AREA>/<SCENARIO>/<cut>/similarity_<pixel size>m (similarity_<pixel size>m_<area>
for another area; see the wiki, Similarity and landscape metrics): sim_<id>_hft<h>.tif, similar_<id>_hft<h>.tif, hri_results.csv, thdi_res.csv, the plots,
similarity_metadata.txt and similarity_run.json (written last).
"""

import argparse
import ast
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime

import numpy as np
from osgeo import gdal, ogr, osr

import hft_state

gdal.UseExceptions()
ogr.UseExceptions()
osr.UseExceptions()

HERE = hft_state.HERE
ROOT = os.path.dirname(HERE)
INVARS = hft_state.INVARS
RESULTS = os.path.join(HERE, 'results')
ECO_FOLDER = 'ecoregion'
ECO_POLYGONS = 'ecoregions.gpkg'
STUDY_AT_AREA = 'study_area'  # subfolder of an area with the study area at its pixel (get_vars_gee.py)
THDI_R = os.path.join(ROOT, 'thdi.R')
PLOTS_R = os.path.join(HERE, 'similarity_plots.R')
RUN_JSON = 'similarity_run.json'
RUN_TXT = 'similarity_metadata.txt'
HRI_CSV = 'hri_results.csv'
THDI_CSV = 'thdi_res.csv'
RUN_DIR_RE = re.compile(r'^similarity_\d+(?:_\d+)?m(?:_[A-Za-z0-9_-]+)?$')  # as run() names them
SHP_RE = re.compile(r'^park_segm_(?P<id>-?\d+)_(?P<k>\d+)(?:_(?P<cut>p\d+))?(?P<optim>_optim)?_class\.shp$')
CUT_RE = re.compile(r'^p\d+$')
HFT_FIELD = 'hclst_m'

NODATA = -9999.0  # of the similarity maps
NO_ECO = -2147483647  # pixels outside the polygons of the ecoregions
CAT_NODATA = 255  # of the maps of the similar areas
P_MIN = 0.001  # similarities (probabilities) of P_MIN or less are 0, as in eHabitat+
MIN_VALID_FRACTION = 0.5  # a coarser pixel needs data in this fraction of its area
RIDGES = (0.0, 1e-8, 1e-6, 1e-4, 1e-2)  # added to the correlation matrix only if it is not positive definite
CHUNK = 1 << 21  # pixels per block of the Mahalanobis distances
BLOCK_BYTES = 256 * 2 ** 20  # memory of the blocks of rows read from the variables
FACTORS = (1, 2, 5, 10)  # pixel sizes offered by the plan: multiples of the pixel of the ecoregion variables
DEFAULTS = {'ECOSIM_AREA': '', 'ECOSIM_CUT': '', 'ECOSIM_RESOLUTION': 0, 'ECOSIM_STAT': 'mean',
            'ECOSIM_NEIGHBOURS': 8, 'ECOSIM_MIN_COVER': 0.5, 'ECOSIM_ECOREGIONS': '', 'ECOSIM_MARINE': ''}
GTIFF = ['TILED=YES', 'COMPRESS=DEFLATE', 'NUM_THREADS=ALL_CPUS', 'BIGTIFF=IF_SAFER']


# Options and inputs -------------------------------------------------------------------------------------------------

def option(cfg, name):
	"""Option NAME of config.py, or its default when the config.py is older than this script."""
	v = cfg.get(name, DEFAULTS[name])
	return DEFAULTS[name] if v is None else v


def rel(path):
	"""PATH relative to the repository, with '/' (as written in the records of the results); absolute if it is
	outside the repository."""
	if not path:
		return None
	if not _inside(path, ROOT):
		return os.path.abspath(path)
	return os.path.relpath(path, ROOT).replace(os.sep, '/')


def user_path(path):
	"""A path of config.py or of the command line: ~ expanded, and relative to the repository if it is not absolute."""
	path = os.path.expanduser(str(path or '').strip())
	if not path:
		return ''
	return path if os.path.isabs(path) else os.path.normpath(os.path.join(ROOT, path))


def variable_key(path):
	"""Name of a variable to pair the two scales: the file name without its order prefix (1_, 2_...) nor extension,
	in lower case (as app/ehabitat_app/inputs.py)."""
	return re.sub(r'^\d+_', '', os.path.splitext(os.path.basename(path))[0]).lower()


def fmt_res(res):
	return f'{res:.0f}' if abs(res - round(res)) < 1e-6 else f'{res:g}'.replace('.', '_')


def scenario_root(study_area, scenario):
	root = os.path.join(RESULTS, str(study_area), str(scenario))
	if not os.path.isdir(root):
		raise ValueError(f'no results of the HFT model in {rel(root)}: run run_hft.py first (or choose another '
		                 f'scenario)')
	return root


def scenario_variables(root, cfg):
	"""ENV_VARS of the scenario (processing_metadata.txt of moranvar_plots.R), else those of config.py."""
	meta = os.path.join(root, 'processing_metadata.txt')
	if os.path.isfile(meta):
		text = open(meta, encoding='utf-8', errors='replace').read()
		m = re.search(r'^ENV_VARS = \{(.*?)^\}', text, re.M | re.S)
		if m:
			try:
				return ast.literal_eval('{' + m.group(1) + '}'), 'processing_metadata.txt'
			except (ValueError, SyntaxError):
				pass
	return dict(cfg.get('ENV_VARS') or {}), 'config.py'


def files_by_key(folder):
	"""{variable key: GeoTIFF} of FOLDER (not of its subfolders)."""
	out = {}
	if os.path.isdir(folder):
		for f in sorted(os.listdir(folder)):
			if f.lower().endswith(hft_state.RASTER_EXT) and not f.startswith('.'):
				out.setdefault(variable_key(f), os.path.join(folder, f))
	return out


def variable_pairs(env_vars, area_dir):
	"""Every variable of the HFTs with its file in the folder AREA_DIR of the ecoregion or comparison area (same name
	without the number in front), in the key 'ecoregion'."""
	files = files_by_key(area_dir)
	out = []
	for name, spec in env_vars.items():
		local = os.path.normpath(os.path.join(INVARS, str(spec.get('file', ''))))
		key = variable_key(local)
		out.append({'name': name, 'key': key, 'local': local, 'ecoregion': files.get(key)})
	return out


def default_area(env_vars):
	"""The folder ecoregion next to the variables of the HFTs."""
	first = next(iter(env_vars.values()), {})
	return os.path.normpath(os.path.join(INVARS, os.path.dirname(str(first.get('file', ''))), ECO_FOLDER))


def area_folder(value, env_vars):
	"""The folder of the ecoregion or comparison area: ECOSIM_AREA (relative to inVars if it is not absolute), or the
	folder ecoregion of the variables of the HFTs."""
	value = os.path.expanduser(str(value or '').strip())
	if not value:
		return default_area(env_vars)
	path = os.path.normpath(value if os.path.isabs(value) else os.path.join(INVARS, value))
	if not os.path.isdir(path):
		raise ValueError(f'the folder of the comparison area (ECOSIM_AREA) does not exist: {path}')
	return path


def area_tag(area_dir, env_vars):
	"""Suffix of the folders of the results for an area other than the ecoregion of the variables of the HFTs."""
	if os.path.normcase(area_dir) == os.path.normcase(default_area(env_vars)):
		return None
	area_dir = os.path.normpath(area_dir)
	if _inside(area_dir, INVARS):
		r = os.path.relpath(area_dir, INVARS).replace(os.sep, '/')
	else:
		r = os.path.basename(area_dir)
		if r == ECO_FOLDER:
			r = os.path.basename(os.path.dirname(area_dir))
	if r.endswith('/' + ECO_FOLDER):
		r = r[:-len(ECO_FOLDER) - 1]
	return re.sub(r'[^A-Za-z0-9_-]', '_', r).strip('_') or 'area'


def _inside(path, folder):
	try:
		return os.path.commonpath([os.path.abspath(path), os.path.abspath(folder)]) == os.path.abspath(folder)
	except ValueError:
		return False


def _record(folder):
	"""gee_vars.json of a folder of get_vars_gee.py ({} if there is none)."""
	try:
		with open(os.path.join(folder, 'gee_vars.json'), encoding='utf-8') as f:
			return json.load(f)
	except (OSError, ValueError):
		return {}


def same_spec(a, b):
	"""Whether two specs of get_vars_gee.py are the same recipe: keys that only one of them has count only if they are
	not null (newer versions of monitoreo add keys such as where: None)."""
	if not isinstance(a, dict) or not isinstance(b, dict):
		return a == b
	return all(a.get(k) == b.get(k) for k in set(a) | set(b) if a.get(k) is not None or b.get(k) is not None)


def area_info(area_dir):
	"""What the folder AREA_DIR holds: RESOLVE ecoregions or a comparison area of the user (get_vars_gee.py --area),
	its name and polygons."""
	rec = _record(area_dir)
	poly = os.path.join(area_dir, ECO_POLYGONS)
	poly = poly if os.path.isfile(poly) else None
	kind, name = 'ecoregions', None
	if rec.get('comparison_area'):
		kind, name = 'comparison area', rec['comparison_area'].get('name')
	elif poly:
		try:
			feats = ecoregion_features(poly)
			ds = ogr.Open(poly)
			lyr = ds.GetLayer(0)
			biomes = {f.GetField('BIOME_NAME') for f in lyr} if lyr.GetLayerDefn().GetFieldIndex('BIOME_NAME') >= 0 else set()
			if 'Comparison area' in biomes or any((x['id'] or 0) < 0 for x in feats['features']):
				kind = 'comparison area'
			name = '; '.join(str(x['name']) for x in feats['features'] if x['name']) or None
		except RuntimeError:
			pass
	elif rec.get('ecoregions'):
		name = '; '.join(str(e.get('ECO_NAME')) for e in rec['ecoregions'] if e.get('ECO_NAME')) or None
	pending = rec.get('comparison_area_pending')  # get_vars_gee.py is replacing the area of the folder: unfinished
	pending = (pending.get('name') or ('the RESOLVE ecoregions' if pending.get('resolve') else 'another area')) \
		if isinstance(pending, dict) else None
	return {'folder': area_dir, 'kind': kind, 'name': name, 'polygons': poly, 'pending': pending}


def reference_folder(area_dir, pairs, grid):
	"""The folder study_area inside AREA_DIR (the study area at the pixel of the area, downloaded by get_vars_gee.py
	when the area does not contain it), if it has every variable on the same lattice as GRID and, when both have
	records of get_vars_gee.py, the same recipe; else None (the HFT statistics come from the area itself). Returns
	(folder or None, files by key, why not)."""
	ref = os.path.join(area_dir, STUDY_AT_AREA)
	if not os.path.isdir(ref):
		return None, None, None
	files = files_by_key(ref)
	missing = [p['key'] for p in pairs if p['key'] not in files]
	if missing:
		return None, None, f"{rel(ref)} has no file of {', '.join(missing)}"
	try:
		g, _ = native_grid([files[p['key']] for p in pairs])
	except (ValueError, RuntimeError) as e:
		return None, None, str(e)
	if not same_lattice(g, grid):
		return None, None, f'{rel(ref)} is not on the grid lines of the area (CRS, pixel or edges)'
	ra, rr = _record(area_dir).get('files') or {}, _record(ref).get('files') or {}
	for p in pairs:
		fa, fr = os.path.basename(p['ecoregion']), os.path.basename(files[p['key']])
		sa, sr = (ra.get(fa) or {}).get('spec'), (rr.get(fr) or {}).get('spec')
		if sa is not None and sr is not None and not same_spec(sa, sr):
			return None, None, f'{fr} of {rel(ref)} was not computed with the same recipe as {fa} of the area'
	return ref, files, None


def find_targets(root, cut, ids=None):
	"""The HFT shapefiles to compare: one per ID, of the cut CUT ('' = the suggested one of each ID; 'results' = the
	scenario folder itself, with SIMIL 1 or 2)."""
	cut = str(cut or '').strip()
	cut_dirs = sorted((d for d in os.listdir(root) if CUT_RE.match(d) and os.path.isdir(os.path.join(root, d))),
	                  key=lambda d: int(d[1:]))
	dirs = [(c, os.path.join(root, c)) for c in cut_dirs] or [('results', root)]
	if cut and cut not in [c for c, _ in dirs]:
		raise ValueError(f"no cut '{cut}' in {rel(root)} (cuts: {', '.join(c for c, _ in dirs)})")
	found = {}
	for c, d in dirs:
		for f in sorted(os.listdir(d)):
			m = SHP_RE.match(f)
			if m:
				found.setdefault(m.group('id'), []).append({
					'id': m.group('id'), 'cut': c, 'k': m.group('k'), 'suggested': bool(m.group('optim')),
					'shp': os.path.join(d, f), 'dir': d,
					'colours': _existing(os.path.join(d, 'hft_colours_' + f[len('park_segm_'):-len('_class.shp')] + '.csv'))})
	if not found:
		raise ValueError(f'no HFT polygons (park_segm_*_class.shp) in {rel(root)}')
	want = [str(i) for i in ids] if ids else sorted(found, key=int)
	out = []
	for i in want:
		if i not in found:
			raise ValueError(f'no HFTs of the ID {i} in {rel(root)}')
		opts = found[i]
		if cut:
			opts = [t for t in opts if t['cut'] == cut]
		elif len(opts) > 1:
			opts = [t for t in opts if t['suggested']]
			if not opts:
				raise ValueError(f'the ID {i} has no suggested cut: choose a cut (ECOSIM_CUT or --cut)')
		if not opts:
			raise ValueError(f'the ID {i} has no HFTs in the cut {cut}')
		out.append(opts[0])
	return out


def _existing(path):
	return path if os.path.isfile(path) else None


def ecoregions_file(cfg_value, area, env_vars):
	"""The polygons of the ecoregions (for the ecoregion of each HFT and TerrEco of the THDI): ECOSIM_ECOREGIONS, or
	those of the area if it holds ecoregions, or those of the ecoregion folder of the variables of the HFTs."""
	if str(cfg_value or '').strip():
		f = user_path(cfg_value)
		if not os.path.isfile(f):
			raise ValueError(f'the file of the ecoregions does not exist: {f}')
		return f
	if area['kind'] == 'ecoregions' and area['polygons']:
		return area['polygons']
	own = area_info(default_area(env_vars))  # those of the ecoregion grid of the variables of the HFTs
	return own['polygons'] if own['kind'] == 'ecoregions' else None


def ecoregion_features(path):
	"""ID field and features (id, name) of a file of ecoregions: ECO_ID, else its first integer field, else the FID."""
	ds = ogr.Open(path)
	lyr = ds.GetLayer(0)
	defn = lyr.GetLayerDefn()
	names = [defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())]
	field = 'ECO_ID' if 'ECO_ID' in names else next(
		(defn.GetFieldDefn(i).GetName() for i in range(defn.GetFieldCount())
		 if defn.GetFieldDefn(i).GetType() in (ogr.OFTInteger, ogr.OFTInteger64)), None)
	label = next((n for n in ('ECO_NAME', 'NAME', 'name', 'nombre') if n in names), None)
	feats = []
	for f in lyr:
		fid = f.GetField(field) if field else f.GetFID() + 1  # + 1: the FIDs of shapefiles start at 0 (0 = none)
		feats.append({'id': int(fid) if fid is not None else None, 'name': f.GetField(label) if label else None})
	return {'file': path, 'field': field or 'FID', 'features': feats}


# Grids -----------------------------------------------------------------------------------------------------------------

class Grid:
	"""A north-up grid: geotransform, size and CRS (WKT)."""

	def __init__(self, gt, width, height, wkt):
		self.gt, self.width, self.height, self.wkt = tuple(float(v) for v in gt), int(width), int(height), wkt
		self.res = abs(self.gt[1])

	@property
	def extent(self):
		x0, y0 = self.gt[0], self.gt[3]
		return [x0, y0 - self.height * self.res, x0 + self.width * self.res, y0]

	def srs(self):
		s = osr.SpatialReference()
		s.ImportFromWkt(self.wkt)
		s.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		return s

	def crs_label(self):
		s = self.srs()
		code = s.GetAuthorityCode(None)
		return f'{s.GetAuthorityName(None)}:{code}' if code else (s.GetName() or 'unknown CRS')


def native_grid(files):
	"""Grid and NoData of the GeoTIFFs FILES (one per variable, all of them on the same grid)."""
	grid, nodata = None, []
	for path in files:
		ds = gdal.Open(path)
		gt = ds.GetGeoTransform()
		if abs(gt[2]) > 1e-9 or abs(gt[4]) > 1e-9 or gt[5] >= 0:
			raise ValueError(f'{os.path.basename(path)} is not a north-up grid')
		if abs(abs(gt[1]) - abs(gt[5])) > 1e-6 * abs(gt[1]):
			raise ValueError(f'the pixels of {os.path.basename(path)} are not square')
		g = Grid(gt, ds.RasterXSize, ds.RasterYSize, ds.GetProjection())
		if grid is None:
			grid = g
			if not g.wkt or g.srs().IsGeographic():
				raise ValueError(f'the variables of {rel(os.path.dirname(path))} need a projected CRS (metres)')
		elif (g.width, g.height) != (grid.width, grid.height) or any(abs(a - b) > 1e-6 * grid.res for a, b in
		                                                          zip(g.gt, grid.gt)) or not g.srs().IsSame(grid.srs()):
			raise ValueError(f'the variables of {rel(os.path.dirname(path))} are not all on the same grid: '
			                 f'{os.path.basename(path)} differs from {os.path.basename(files[0])}')
		nodata.append(ds.GetRasterBand(1).GetNoDataValue())
	return grid, nodata


def same_lattice(a, b):
	"""Whether the grids A and B have the same CRS and pixel, with their edges on the same lines."""
	if abs(a.res - b.res) > 1e-6 * a.res or not a.srs().IsSame(b.srs()):
		return False
	return all(abs(d / a.res - round(d / a.res)) < 1e-6 for d in (a.gt[0] - b.gt[0], a.gt[3] - b.gt[3]))


class CoarseGrid(Grid):
	"""The grid of the similarity: blocks of FACTOR x FACTOR pixels of the native grid, with edges at multiples of
	the pixel size when the native grid has its edges at multiples of its own pixel (as the grids of MonitorEO), so
	that every native pixel is inside one pixel of the coarse grid."""

	def __init__(self, native, res):
		f = res / native.res
		if res <= 0 or f < 1 - 1e-6 or abs(f - round(f)) > 1e-6:
			raise ValueError(f'the pixel size of the similarity must be a multiple of the pixel of the ecoregion '
			                 f'variables ({fmt_res(native.res)} m): {fmt_res(res)} m is not')
		f = int(round(f))
		x0, y0 = native.gt[0], native.gt[3]
		offx = offy = 0
		if f > 1 and all(abs(v / native.res - round(v / native.res)) < 1e-6 for v in (x0, y0)):
			cx0 = math.floor(x0 / res + 1e-9) * res
			cy0 = math.ceil(y0 / res - 1e-9) * res
			offx, offy = int(round((x0 - cx0) / native.res)), int(round((cy0 - y0) / native.res))
		else:
			cx0, cy0 = x0, y0
		width = -(-(offx + native.width) // f)
		height = -(-(offy + native.height) // f)
		super().__init__((cx0, native.res * f, 0, cy0, 0, -native.res * f), width, height, native.wkt)
		self.native, self.factor, self.offx, self.offy = native, f, offx, offy


def read_coarse(bands, nodata, cg, r0, c0, nr, nc):
	"""Values of every band on the window (R0, C0, NR, NC) of the coarse grid CG: the native values or, with a coarser
	pixel, the mean of its native pixels with data if they cover at least MIN_VALID_FRACTION of it. NaN = no data."""
	f, n = cg.factor, cg.native
	# native rows and columns of the window (the coarse grid may start before the native one: offx, offy)
	nr0, nc0 = r0 * f - cg.offy, c0 * f - cg.offx
	rr0, rr1 = max(0, nr0), min(n.height, nr0 + nr * f)
	cc0, cc1 = max(0, nc0), min(n.width, nc0 + nc * f)
	out = []
	for band, nd in zip(bands, nodata):
		buf = np.full((nr * f, nc * f), np.nan, dtype=np.float32)
		if rr1 > rr0 and cc1 > cc0:
			a = band.ReadAsArray(cc0, rr0, cc1 - cc0, rr1 - rr0).astype(np.float32, copy=False)
			bad = ~np.isfinite(a)
			if nd is not None and np.isfinite(nd):
				bad |= a == np.float32(nd)
			a[bad] = np.nan
			buf[rr0 - nr0:rr1 - nr0, cc0 - nc0:cc1 - nc0] = a
		if f == 1:
			out.append(buf)
			continue
		b = buf.reshape(nr, f, nc, f)
		cnt = np.isfinite(b).sum(axis=(1, 3))
		with np.errstate(invalid='ignore', divide='ignore'):
			mean = np.nansum(b, axis=(1, 3), dtype=np.float64) / cnt
		mean[cnt < MIN_VALID_FRACTION * f * f] = np.nan
		out.append(mean.astype(np.float32))
	return out


def open_bands(files):
	dss = [gdal.Open(f) for f in files]
	return dss, [d.GetRasterBand(1) for d in dss]


def load_stack(files, nodata, cg):
	"""The pixels of the coarse grid with data in every variable (complete cases): their mask (rows x columns) and
	their values (pixels x variables, float32), read in blocks of rows."""
	dss, bands = open_bands(files)
	p = len(bands)
	rows = max(1, int(BLOCK_BYTES // (cg.width * cg.factor * cg.factor * 4 * (p + 1))))
	valid = np.zeros((cg.height, cg.width), dtype=bool)
	# two passes (the mask, then the values) so that the values are not held twice in memory
	for r0 in range(0, cg.height, rows):
		nr = min(rows, cg.height - r0)
		ok = np.ones((nr, cg.width), dtype=bool)
		for v in read_coarse(bands, nodata, cg, r0, 0, nr, cg.width):
			ok &= np.isfinite(v)
		valid[r0:r0 + nr] = ok
	X = np.empty((int(valid.sum()), p), dtype=np.float32)
	i = 0
	for r0 in range(0, cg.height, rows):
		nr = min(rows, cg.height - r0)
		ok = valid[r0:r0 + nr]
		n = int(ok.sum())
		if n:
			for j, v in enumerate(read_coarse(bands, nodata, cg, r0, 0, nr, cg.width)):
				X[i:i + n, j] = v[ok]
			i += n
	return valid, X


# HFTs on the grid of the similarity ----------------------------------------------------------------------------------

def _mem_vector():
	"""An in-memory vector dataset: driver MEM since GDAL 3.11 (Memory before, deprecated since then)."""
	drv = gdal.GetDriverByName('MEM')
	if drv is not None and drv.GetMetadataItem(gdal.DCAP_VECTOR) == 'YES':
		return drv.Create('', 0, 0, 0, gdal.GDT_Unknown)
	return ogr.GetDriverByName('Memory').CreateDataSource('')


def _layer_in(path, srs, extent=None):
	"""The polygons of PATH in the CRS SRS (an in-memory copy, with the same feature IDs, if they have another one),
	with their datasource. EXTENT (xmin, ymin, xmax, ymax in SRS): only the polygons near it, so that a global file
	(e.g. all the RESOLVE ecoregions) is not reprojected far from the grid, where its CRS is not defined."""
	ds = ogr.Open(path)
	lyr = ds.GetLayer(0)
	src = lyr.GetSpatialRef()
	if src is None:
		return ds, lyr
	src = src.Clone()
	src.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	if extent is not None:
		w = max(extent[2] - extent[0], extent[3] - extent[1]) * 0.05
		box = (extent[0] - w, extent[1] - w, extent[2] + w, extent[3] + w)
		if not src.IsSame(srs):
			box = osr.CoordinateTransformation(srs, src).TransformBounds(*box, 21)
		lyr.SetSpatialFilterRect(*box)
	if src.IsSame(srs):
		return ds, lyr
	tr = osr.CoordinateTransformation(src, srs)
	mem = _mem_vector()
	out = mem.CreateLayer('t', srs=srs, geom_type=ogr.wkbMultiPolygon)
	defn = lyr.GetLayerDefn()
	for i in range(defn.GetFieldCount()):
		out.CreateField(defn.GetFieldDefn(i))
	for f in lyr:
		g = f.GetGeometryRef()
		if g is None:
			continue
		g = g.Clone()
		try:
			g.Transform(tr)
		except RuntimeError:  # outside the domain of the CRS of the grid: it cannot overlap it
			continue
		nf = ogr.Feature(out.GetLayerDefn())
		nf.SetFrom(f)
		nf.SetFID(f.GetFID())
		nf.SetGeometry(g)
		out.CreateFeature(nf)
	return (mem, ds), out


def hft_areas(shp):
	"""Area (ha) of the polygons of every HFT, in the CRS of the shapefile (or in metres of its UTM if geographic)."""
	ds = ogr.Open(shp)
	lyr = ds.GetLayer(0)
	srs = lyr.GetSpatialRef()
	tr = None
	if srs is not None and srs.IsGeographic():
		srs = srs.Clone()
		srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		x0, x1, y0, y1 = lyr.GetExtent()
		utm = osr.SpatialReference()
		utm.SetWellKnownGeogCS('WGS84')
		utm.SetUTM(int(((x0 + x1) / 2 + 180) // 6) + 1, (y0 + y1) / 2 >= 0)
		utm.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
		tr = osr.CoordinateTransformation(srs, utm)
	unit = srs.GetLinearUnits() if srs is not None and not srs.IsGeographic() else 1.0
	areas = {}
	for f in lyr:
		g = f.GetGeometryRef()
		h = f.GetField(HFT_FIELD)
		if g is None or h is None:
			continue
		if tr:
			g = g.Clone()
			g.Transform(tr)
		areas[int(h)] = areas.get(int(h), 0.0) + g.GetArea() * unit * unit / 1e4
	return areas


def hft_cover(shp, cg, local_res, min_cover):
	"""The HFTs of SHP on the coarse grid: in a window around them, the HFT of every pixel (the one that covers most
	of it, if it covers at least MIN_COVER of it; 0 = none) and the fraction it covers, from the polygons drawn on a
	grid SUB times finer (SUB = pixel of the similarity / pixel of the HFT map, so that with nested grids the
	fractions are those of the HFT map)."""
	src, lyr = _layer_in(shp, cg.srs())
	if lyr.GetLayerDefn().GetFieldIndex(HFT_FIELD) < 0:
		raise ValueError(f'{os.path.basename(shp)} has no column of the HFTs ({HFT_FIELD})')
	hfts = sorted({int(f.GetField(HFT_FIELD)) for f in lyr if f.GetField(HFT_FIELD) is not None})
	lyr.ResetReading()
	xmin, xmax, ymin, ymax = lyr.GetExtent()
	x0, y0, res = cg.gt[0], cg.gt[3], cg.res
	gx0, gy0, gx1, gy1 = cg.extent
	clipped = xmin < gx0 - 1e-6 or xmax > gx1 + 1e-6 or ymin < gy0 - 1e-6 or ymax > gy1 + 1e-6
	c0 = max(0, int(math.floor((xmin - x0) / res)) - 1)
	c1 = min(cg.width, int(math.ceil((xmax - x0) / res)) + 1)
	r0 = max(0, int(math.floor((y0 - ymax) / res)) - 1)
	r1 = min(cg.height, int(math.ceil((y0 - ymin) / res)) + 1)
	if c1 <= c0 or r1 <= r0:
		raise ValueError(f'the HFTs of {os.path.basename(shp)} are outside the grid of the ecoregion variables')
	nr, nc = r1 - r0, c1 - c0
	sub = max(1, int(math.ceil(res / local_res - 1e-6))) if local_res else 10
	while sub > 1 and nr * nc * sub * sub > 64e6:
		sub -= 1
	mem = gdal.GetDriverByName('MEM').Create('', nc * sub, nr * sub, 1, gdal.GDT_Int32)
	mem.SetGeoTransform((x0 + c0 * res, res / sub, 0, y0 - r0 * res, 0, -res / sub))
	mem.SetProjection(cg.wkt)
	gdal.RasterizeLayer(mem, [1], lyr, options=[f'ATTRIBUTE={HFT_FIELD}'])
	a = mem.GetRasterBand(1).ReadAsArray().reshape(nr, sub, nc, sub).transpose(0, 2, 1, 3).reshape(nr * nc, sub * sub)
	counts = np.stack([(a == h).sum(axis=1) for h in hfts], axis=1) if hfts else np.zeros((nr * nc, 0))
	frac = counts / float(sub * sub)
	best = frac.argmax(axis=1) if hfts else np.zeros(nr * nc, dtype=int)
	cover = frac[np.arange(nr * nc), best] if hfts else np.zeros(nr * nc)
	assign = np.where((cover >= min_cover - 1e-9) & (cover > 0), np.asarray(hfts)[best] if hfts else 0, 0)
	return {'hfts': hfts, 'window': (r0, c0, nr, nc), 'assign': assign.reshape(nr, nc).astype(np.int32),
	        'cover': cover.reshape(nr, nc), 'sub': sub, 'clipped': clipped}


def local_resolution(pairs):
	"""Pixel of the variables of the HFT map (the local scale)."""
	for p in pairs:
		if os.path.isfile(p['local']):
			return abs(gdal.Open(p['local']).GetGeoTransform()[1])
	return None


def ecoregion_ids(eco, cg, valid):
	"""ECO_ID (or the ID field) of every pixel with data, when the file has several ecoregions (else None)."""
	if not eco or len(eco['features']) < 2:
		return None
	src, lyr = _layer_in(eco['file'], cg.srs(), cg.extent)
	mem = gdal.GetDriverByName('MEM').Create('', cg.width, cg.height, 1, gdal.GDT_Int32)
	mem.SetGeoTransform(cg.gt)
	mem.SetProjection(cg.wkt)
	mem.GetRasterBand(1).Fill(NO_ECO)  # not 0: ECO_ID 0 is an ecoregion of RESOLVE (Rock and Ice)
	if eco['field'] != 'FID':
		gdal.RasterizeLayer(mem, [1], lyr, options=[f'ATTRIBUTE={eco["field"]}'])
	else:  # without an ID field, the ID of every ecoregion is its feature ID + 1 (as in ecoregion_features)
		for f in lyr:
			g = f.GetGeometryRef()
			if g is None:
				continue
			tmp = _mem_vector()
			t = tmp.CreateLayer('t', srs=cg.srs(), geom_type=ogr.wkbMultiPolygon)
			nf = ogr.Feature(t.GetLayerDefn())
			nf.SetGeometry(g.Clone())
			t.CreateFeature(nf)
			gdal.RasterizeLayer(mem, [1], t, burn_values=[f.GetFID() + 1])
	return mem.GetRasterBand(1).ReadAsArray()[valid]


# Similarity -----------------------------------------------------------------------------------------------------------

def _compact_index(valid, rows, cols):
	"""Positions in the list of pixels with data (load_stack) of the pixels (ROWS, COLS), which must have data."""
	row_start = np.concatenate([[0], np.cumsum(valid.sum(axis=1))])
	out = np.empty(len(rows), dtype=np.int64)
	for r in np.unique(rows):
		sel = rows == r
		cs = np.cumsum(valid[r])
		out[sel] = row_start[r] + cs[cols[sel]] - 1
	return out


def reference(Xh, eco_sd):
	"""Mean, standard deviations and Cholesky factor of the correlation matrix of the values XH (pixels x variables)
	of the HFT, without the variables that are constant in it (relative to their deviation in the ecoregion); RIDGE is
	what had to be added to the diagonal of the correlation matrix to make it positive definite (0 almost always)."""
	from scipy.linalg import cholesky
	Xh = np.asarray(Xh, dtype=np.float64)
	sd = Xh.std(axis=0, ddof=1)
	keep = np.flatnonzero(sd > 1e-6 * np.maximum(eco_sd, 1e-30))
	if keep.size == 0:
		return None
	mu, sd_k = Xh[:, keep].mean(axis=0), sd[keep]
	corr = np.atleast_2d(np.corrcoef((Xh[:, keep] - mu) / sd_k, rowvar=False))
	for ridge in RIDGES:
		try:
			L = cholesky(corr + ridge * np.eye(len(keep)), lower=True, check_finite=False)
			return {'keep': keep, 'mu': mu, 'sd': sd_k, 'L': L, 'ridge': ridge}
		except np.linalg.LinAlgError:
			continue
	return None


def similarity_values(X, ref, domain=None):
	"""Similarity of the pixels DOMAIN (indices in X; None = all) with the HFT: chi-square probability (degrees of
	freedom = the variables used) of the squared Mahalanobis distance, 0 where it is P_MIN or less."""
	from scipy.linalg import solve_triangular
	from scipy.special import chdtrc
	keep, mu, sd, L = ref['keep'], ref['mu'], ref['sd'], ref['L']
	n = len(X) if domain is None else len(domain)
	out = np.empty(n, dtype=np.float32)
	for s in range(0, n, CHUNK):
		rows = X[s:s + CHUNK] if domain is None else X[domain[s:s + CHUNK]]
		z = (rows[:, keep].astype(np.float64) - mu) / sd
		w = solve_triangular(L, z.T, lower=True, check_finite=False)
		d2 = np.einsum('ij,ij->j', w, w)
		out[s:s + CHUNK] = chdtrc(len(keep), d2)
	out[out <= P_MIN] = 0
	return out


def write_raster(path, arr, cg, dtype, nodata, resampling):
	"""ARR (rows x columns of CG) as a compressed, tiled GeoTIFF with internal overviews."""
	drv = gdal.GetDriverByName('GTiff')
	pred = 'PREDICTOR=3' if dtype == gdal.GDT_Float32 else 'PREDICTOR=2'
	ds = drv.Create(path, cg.width, cg.height, 1, dtype, options=GTIFF + [pred])
	ds.SetGeoTransform(cg.gt)
	ds.SetProjection(cg.wkt)
	b = ds.GetRasterBand(1)
	b.SetNoDataValue(nodata)
	b.WriteArray(arr)
	levels, n = [], 2
	while max(cg.width, cg.height) / n >= 256:
		levels.append(n)
		n *= 2
	if levels:
		ds.BuildOverviews(resampling, levels)
	ds.FlushCache()
	ds = None


def landscape(similar, structure, n_hft):
	"""Patches of the similar pixels (boolean rows x columns): their labels and the size (pixels) of every patch."""
	from scipy import ndimage
	labels, n = ndimage.label(similar, structure=structure)
	sizes = np.bincount(labels.ravel(), minlength=n + 1)[1:] if n else np.zeros(0, dtype=np.int64)
	return labels, sizes


def run_hft(h, ctx, target, Xh, cells, eco_v, eco_of_hft, outdir):
	"""Similarity, maps and landscape metrics of the HFT H of TARGET, with XH the values of its pixels (pixels x
	variables, on the grid of the similarity, from the area or from the study area at its pixel) and CELLS the indices
	(row * width + column) of its pixels in the grid of the area, if it is there (left out of the similar areas).
	Returns its row of hri_results.csv and its record."""
	t0 = time.time()
	X, valid, cg = ctx['X'], ctx['valid'], ctx['cg']
	n = len(Xh)
	hid = target['id']
	rec = {'id': hid, 'hft': int(h), 'ecoregion': eco_of_hft, 'status': 'ok', 'ref_pixels': int(n), 'vars_used': None,
	       'vars_dropped': [], 'ridge': None, 'area_ha': round(ctx['areas'][h], 2) if h in ctx['areas'] else None,
	       'threshold': None, 'sim': None, 'similar': None}
	row = {'ecoregion': eco_of_hft, 'wdpaid': hid, 'wdpaid_hft': f'{hid}_{h}', 'hft': int(h),
	       'HFTArea_ha': rec['area_ha'], 'PixelSize_m': cg.res, 'NrPixHFT': n}
	for j, pr in enumerate(ctx['pairs']):
		v = Xh[:, j] if n else np.zeros(0)
		row[f"{pr['key']}_pamin"] = float(v.min()) if n else None
		row[f"{pr['key']}_pamax"] = float(v.max()) if n else None
		row[f"{pr['key']}_pamean"] = float(v.mean()) if n else None
	if n < ctx['min_pixels']:
		rec['status'] = row['Status'] = 'too_few_pixels'
		print(f'  HFT {h}: {n} pixels of {fmt_res(cg.res)} m (at least {ctx["min_pixels"]} are needed): not computed',
		      flush=True)
		return row, rec
	ref = reference(Xh, ctx['eco_sd'])
	if ref is None:  # every variable is constant in the HFT: no distance can be computed
		rec['status'] = row['Status'] = 'too_few_pixels'
		rec['vars_dropped'] = [pr['key'] for pr in ctx['pairs']]
		print(f'  HFT {h}: all its variables are constant in its {n} pixels: not computed', flush=True)
		return row, rec
	rec['vars_used'] = row['NrVars'] = int(len(ref['keep']))
	rec['ridge'] = ref['ridge']
	# variables constant in the HFT: their covariance would be singular (the failure that the random noise of
	# ehab_optim.py avoided), so they are left out of the distance (one degree of freedom less each)
	rec['vars_dropped'] = [pr['key'] for j, pr in enumerate(ctx['pairs']) if j not in set(ref['keep'].tolist())]
	# pixels compared: those with data in the area (in the ecoregion of the HFT when the area has several ecoregions)
	domain = np.flatnonzero(eco_v == eco_of_hft) if eco_v is not None and eco_of_hft is not None else None
	sim_v = similarity_values(X, ref, domain)
	sim_full = np.full(cg.height * cg.width, NODATA, dtype=np.float32)
	if domain is None:
		sim_full[valid.ravel()] = sim_v
	else:
		sim_full[np.flatnonzero(valid.ravel())[domain]] = sim_v
	del sim_v, domain
	s_ref = similarity_values(Xh, ref)  # the similarity of the pixels of the HFT itself
	pos = s_ref[s_ref > 0]
	status = 'few_pixels' if n < ctx['few_pixels'] else 'ok'
	agg = np.median if ctx['stat'] == 'median' else np.mean
	sim_name = f'sim_{hid}_hft{h}.tif'
	write_raster(os.path.join(outdir, sim_name), sim_full.reshape(cg.height, cg.width), cg, gdal.GDT_Float32,
	             NODATA, 'AVERAGE')
	rec['sim'] = sim_name
	cat = np.full(cg.height * cg.width, CAT_NODATA, dtype=np.uint8)
	cat[sim_full != NODATA] = 0
	if not pos.size:
		status = 'no_reference_similarity'
		row.update({'Status': status, 'NrSimPixEco': None})
	else:
		T = float(agg(pos))
		rec['threshold'] = T
		similar = sim_full >= T
		similar[cells] = False
		del sim_full
		labels, sizes = landscape(similar.reshape(cg.height, cg.width), ctx['structure'], n)
		n_sim = int(similar.sum())
		del similar
		row.update({ctx['sim_label']: T, 'NrSimPixEco': n_sim, 'NrSimPixHFT': int((s_ref >= T).sum()),
		            'NrSimPixEco_NrPixHFT_ratio': n_sim / n, 'NrSimLcpPatEco': int(len(sizes)),
		            'SimAreaEco_km2': n_sim * cg.res * cg.res / 1e6})
		if len(sizes):
			big = sizes >= n
			largest = int(sizes.max())
			hsr = float(sizes[big].sum()) / n
			row.update({'NrPixLargLcpPat_NrPixHFT_ratio': largest / n, 'HSR': hsr if hsr > 0 else largest / n,
			            'NrSimLcpPatEqLargArea': int(big.sum()), 'NrPixLargLcpPat': largest,
			            'Aggregation': n_sim / len(sizes), 'LargPatArea_km2': largest * cg.res * cg.res / 1e6})
			lab = labels.ravel()
			is_big = np.concatenate([[False], big])
			on = lab > 0
			cat[on] = np.where(is_big[lab[on]], 2, 1)
			del on, lab
		else:
			status = 'no_similar'
		del labels
		row['Status'] = status
	cat[cells] = 3
	rec['status'] = status
	sim_cat = f'similar_{hid}_hft{h}.tif'
	write_raster(os.path.join(outdir, sim_cat), cat.reshape(cg.height, cg.width), cg, gdal.GDT_Byte, CAT_NODATA,
	             'NEAREST')
	rec['similar'] = sim_cat
	msg = (f"threshold {rec['threshold']:.3f}, {row.get('NrSimPixEco') or 0} similar pixels in "
	       f"{row.get('NrSimLcpPatEco') or 0} patches, HSR {_num(row.get('HSR'))}" if rec['threshold'] is not None
	       else 'no pixel of the HFT has a similarity above 0')
	ridge = f", ridge {ref['ridge']:g}" if ref['ridge'] else ''
	dropped = f" (left out, constant in the HFT: {', '.join(rec['vars_dropped'])})" if rec['vars_dropped'] else ''
	print(f"  HFT {h}: {n} pixels, {rec['vars_used']} variables{dropped}{ridge}: {msg} ({status}, "
	      f"{time.time() - t0:.0f} s)", flush=True)
	return row, rec


def _num(v):
	return 'NA' if v is None or (isinstance(v, float) and not math.isfinite(v)) else f'{v:.6g}'


HRI_COLUMNS = ['ecoregion', 'wdpaid', 'wdpaid_hft', '{sim}', 'NrSimPixEco', 'NrPixHFT', 'NrSimPixHFT',
               'NrSimPixEco_NrPixHFT_ratio', 'NrSimLcpPatEco', 'NrPixLargLcpPat_NrPixHFT_ratio', 'HSR',
               'NrSimLcpPatEqLargArea', 'NrPixLargLcpPat', 'Aggregation', 'hft', 'HFTArea_ha', 'SimAreaEco_km2',
               'LargPatArea_km2', 'PixelSize_m', 'NrVars', 'Status']


def write_hri(path, rows, pairs, sim_label):
	cols = [c.replace('{sim}', sim_label) for c in HRI_COLUMNS]
	for pr in pairs:
		cols += [f"{pr['key']}_pamin", f"{pr['key']}_pamax", f"{pr['key']}_pamean"]
	with open(path, 'w', encoding='utf-8') as f:
		f.write(' '.join(cols) + '\n')
		for r in rows:
			vals = []
			for c in cols:
				v = r.get(c)
				if v is None or (isinstance(v, float) and not math.isfinite(v)):
					vals.append('NA')
				elif isinstance(v, float):
					vals.append(f'{v:.6g}')
				else:
					vals.append(str(v).replace(' ', '_'))
			f.write(' '.join(vals) + '\n')


# Plan and run ----------------------------------------------------------------------------------------------------------

def _settings(cfg, study_area=None, scenario=None, cut=None, ids=None, res=None, stat=None, neighbours=None,
              min_cover=None, ecoregions=None, marine=None, area=None):
	s = {'study_area': study_area or cfg.get('STUDY_AREA'),
	     'scenario': scenario or hft_state.scenario_name(cfg),
	     'area': option(cfg, 'ECOSIM_AREA') if area is None else area,
	     'cut': option(cfg, 'ECOSIM_CUT') if cut is None else cut,
	     'ids': ids or None,
	     'res': float(option(cfg, 'ECOSIM_RESOLUTION') if res is None else res),
	     'stat': str(option(cfg, 'ECOSIM_STAT') if stat is None else stat).lower(),
	     'neighbours': int(option(cfg, 'ECOSIM_NEIGHBOURS') if neighbours is None else neighbours),
	     'min_cover': float(option(cfg, 'ECOSIM_MIN_COVER') if min_cover is None else min_cover),
	     'ecoregions': option(cfg, 'ECOSIM_ECOREGIONS') if ecoregions is None else ecoregions,
	     'marine': option(cfg, 'ECOSIM_MARINE') if marine is None else marine}
	if s['stat'] not in ('mean', 'median'):
		raise ValueError(f"ECOSIM_STAT must be 'mean' or 'median', not {s['stat']!r}")
	if s['neighbours'] not in (4, 8):
		raise ValueError(f"ECOSIM_NEIGHBOURS must be 8 or 4, not {s['neighbours']}")
	if not 0 < s['min_cover'] <= 1:
		raise ValueError(f"ECOSIM_MIN_COVER must be more than 0 and at most 1, not {s['min_cover']}")
	if s['res'] < 0:
		raise ValueError('ECOSIM_RESOLUTION cannot be negative')
	return s


def _inputs(s, cfg):
	"""Results, variables, area, its grid and the grid of the HFT statistics of the settings S."""
	root = scenario_root(s['study_area'], s['scenario'])
	env_vars, source = scenario_variables(root, cfg)
	if not env_vars:
		raise ValueError(f'no variables (ENV_VARS) for the scenario {s["scenario"]}')
	area_dir = area_folder(s['area'], env_vars)
	pairs = variable_pairs(env_vars, area_dir)
	area = area_info(area_dir)
	area['tag'] = area_tag(area_dir, env_vars)
	if area.get('pending'):
		raise ValueError(f"the download of {area['pending']} into {rel(area_dir)} did not finish, so its files are a mix of "
		                 f"two areas: run get_vars_gee.py again with the same options to finish it (or with the previous "
		                 f"area to go back)")
	inp = {'root': root, 'env_vars': env_vars, 'source': source, 'pairs': pairs, 'area': area,
	       'missing': [p['key'] for p in pairs if not p['ecoregion']], 'native': None, 'nodata': None,
	       'ref': None, 'ref_files': None, 'ref_note': None}
	if inp['missing']:
		return inp
	inp['native'], inp['nodata'] = native_grid([p['ecoregion'] for p in pairs])
	inp['ref'], files, inp['ref_note'] = reference_folder(area_dir, pairs, inp['native'])
	if inp['ref']:
		inp['ref_files'] = [files[p['key']] for p in pairs]
		inp['ref_native'], inp['ref_nodata'] = native_grid(inp['ref_files'])
	inp['eco_file'] = ecoregions_file(s['ecoregions'], area, env_vars)
	return inp


def _area_record(inp):
	a = inp['area']
	return {'folder': rel(a['folder']), 'tag': a['tag'], 'kind': a['kind'], 'name': a['name'],
	        'polygons': rel(a['polygons']), 'reference': rel(inp['ref'] or a['folder']),
	        'reference_note': inp['ref_note']}


def _other_area(outdir, inp):
	"""The area of the results already in OUTDIR when it is another one (two areas whose folders give the same
	name), else None."""
	old = _read_json(os.path.join(outdir, RUN_JSON))
	f = (old.get('area') or {}).get('folder')
	if f is None and old:
		f = rel(default_area(inp['env_vars']))  # runs made before the comparison areas: the default area
	return f if f and f != rel(inp['area']['folder']) else None


def _run_folder(t, cg, tag):
	return os.path.join(t['dir'], f'similarity_{fmt_res(cg.res)}m' + (f'_{tag}' if tag else ''))


def hft_values(t, inp, cg, cg_ref, min_cover, local_res, valid=None, X=None):
	"""The HFTs of the target T on the grid of their statistics (the area, or the study area at its pixel): {hft:
	values (pixels x variables)}, the cover of the area grid (None if the HFTs are outside it) and the HFTs."""
	if inp['ref']:
		try:
			cov_ref = hft_cover(t['shp'], cg_ref, local_res, min_cover)
		except ValueError:
			cov_ref = None
		if cov_ref is None or cov_ref['clipped']:
			raise ValueError(f"the HFTs of the ID {t['id']} are not all inside {rel(inp['ref'])} (the study area at the "
			                 f"pixel of the area): it was downloaded for another study area or buffer; download the "
			                 f"area again with get_vars_gee.py") from None
		dss, bands = open_bands(inp['ref_files'])
		r0, c0, nr, nc = cov_ref['window']
		vals = read_coarse(bands, inp['ref_nodata'], cg_ref, r0, c0, nr, nc)
		ok = np.ones((nr, nc), dtype=bool)
		for v in vals:
			ok &= np.isfinite(v)
		out = {}
		for h in cov_ref['hfts']:
			m = (cov_ref['assign'] == h) & ok
			out[h] = np.stack([v[m] for v in vals], axis=1)
		try:
			cov = hft_cover(t['shp'], cg, local_res, min_cover)
		except ValueError:  # the area does not reach the HFTs
			cov = None
		return out, cov, cov_ref['hfts']
	try:
		cov = hft_cover(t['shp'], cg, local_res, min_cover)
	except ValueError:
		cov = None
	if cov is None or cov['clipped']:
		why = (f"its folder {STUDY_AT_AREA} cannot be used ({inp['ref_note']})" if inp['ref_note'] else
		       f"it has no folder {STUDY_AT_AREA} with the study area at its pixel")
		raise ValueError(f"the HFTs of the ID {t['id']} are {'outside' if cov is None else 'not all inside'} the grid of "
		                 f"{rel(inp['area']['folder'])}, and {why}: get_vars_gee.py downloads it with an ecoregion or "
		                 f"comparison area that does not contain the study area")
	r0, c0, nr, nc = cov['window']
	out = {}
	if X is None:  # plan: the values of the window only
		dss, bands = open_bands([p['ecoregion'] for p in inp['pairs']])
		vals = read_coarse(bands, inp['nodata'], cg, r0, c0, nr, nc)
		ok = np.ones((nr, nc), dtype=bool)
		for v in vals:
			ok &= np.isfinite(v)
		for h in cov['hfts']:
			m = (cov['assign'] == h) & ok
			out[h] = np.stack([v[m] for v in vals], axis=1)
	else:
		for h in cov['hfts']:
			rr, cc = np.nonzero((cov['assign'] == h) & valid[r0:r0 + nr, c0:c0 + nc])
			out[h] = X[_compact_index(valid, rr + r0, cc + c0)] if len(rr) else np.zeros((0, X.shape[1]), np.float32)
	return out, cov, cov['hfts']


def _cells(cov, h, cg):
	"""Indices (row * width + column) of the pixels of the HFT H in the area grid (none if the area is elsewhere)."""
	if cov is None:
		return np.zeros(0, dtype=np.int64)
	r0, c0, nr, nc = cov['window']
	rr, cc = np.nonzero(cov['assign'] == h)
	return (rr + r0).astype(np.int64) * cg.width + (cc + c0)


def memory_mb(cg, p):
	"""Rough peak memory (MB) of a run on the grid CG with P variables (measured: 3.1 GB at 100 m with Guajares)."""
	return int(round(cg.width * cg.height * (4 * p + 14) / 2 ** 20 + 600))


def comparison_areas(study_area=None, scenario=None, cfg=None):
	"""The folders that can be the ecoregion or comparison area of the HFTs of a scenario: every folder ecoregion of
	inVars (up to 3 levels down) and ECOSIM_AREA, with what they are and the variables of the scenario they lack."""
	cfg = cfg or hft_state.read_config()
	root = scenario_root(study_area or cfg.get('STUDY_AREA'), scenario or hft_state.scenario_name(cfg))
	env_vars, _ = scenario_variables(root, cfg)
	default = default_area(env_vars)
	found = []
	for base, dirs, _files in os.walk(INVARS):
		depth = os.path.relpath(base, INVARS).count(os.sep) + (0 if base == INVARS else 1)
		dirs[:] = [d for d in dirs if not d.startswith(('.', '_')) and d != STUDY_AT_AREA]
		if depth >= 3:
			dirs[:] = []
		if os.path.basename(base) == ECO_FOLDER:
			found.append(base)
	conf = option(cfg, 'ECOSIM_AREA')
	if conf:
		try:
			found.append(area_folder(conf, env_vars))
		except ValueError:
			pass
	out, seen = [], set()
	for d in [default] + sorted(found):
		key = os.path.normcase(os.path.normpath(d))
		if key in seen or not os.path.isdir(d):
			continue
		seen.add(key)
		pairs = variable_pairs(env_vars, d)
		info = area_info(d)
		item = {'folder': os.path.relpath(d, INVARS).replace(os.sep, '/') if _inside(d, INVARS) else d,
		        'pending': info.get('pending'),
		        'default': key == os.path.normcase(default), 'kind': info['kind'], 'name': info['name'],
		        'missing': [p['key'] for p in pairs if not p['ecoregion']], 'res': None, 'crs': None, 'width': None,
		        'height': None, 'reference': None, 'polygons': rel(info['polygons'])}
		if not item['missing']:
			try:
				g, _ = native_grid([p['ecoregion'] for p in pairs])
				item.update({'res': g.res, 'crs': g.crs_label(), 'width': g.width, 'height': g.height})
				ref, _, _ = reference_folder(d, pairs, g)
				item['reference'] = rel(ref)
			except (ValueError, RuntimeError) as e:
				item['error'] = str(e)
		out.append(item)
	return out


def plan(study_area=None, scenario=None, cut=None, ids=None, res=None, min_cover=None, ecoregions=None, cfg=None,
         area=None):
	"""What a run would do, without computing it: variables of both scales, grids, the area and the pixels of every
	HFT on the grid of the similarity (see the app, page Similarity). Raises ValueError for the errors of the user."""
	cfg = cfg or hft_state.read_config()
	s = _settings(cfg, study_area, scenario, cut, ids, res, min_cover=min_cover, ecoregions=ecoregions, area=area)
	inp = _inputs(s, cfg)
	pairs = inp['pairs']
	out = {'study_area': s['study_area'], 'scenario': s['scenario'], 'variables_from': inp['source'],
	       'variables': [{k: (rel(v) if k in ('local', 'ecoregion') else v) for k, v in p.items()} for p in pairs],
	       'missing': inp['missing'], 'errors': [], 'warnings': [], 'area': None,
	       'min_pixels': max(5, len(pairs) + 1), 'few_pixels': 3 * len(pairs), 'targets': [], 'ecoregions': None,
	       'native': None, 'res': None, 'factor': None, 'grid': None, 'memory_mb': None, 'resolutions': []}
	if inp['missing']:
		out['errors'].append(f"variables without their file in {rel(inp['area']['folder'])}: "
		                     + ', '.join(inp['missing']))
		return out
	native = inp['native']
	out['area'] = _area_record(inp)
	if inp['ref_note']:
		out['warnings'].append(f"{rel(os.path.join(inp['area']['folder'], STUDY_AT_AREA))} is not used: "
		                       f"{inp['ref_note']}")
	out['native'] = {'res': native.res, 'crs': native.crs_label(), 'width': native.width, 'height': native.height,
	                 'extent': native.extent}
	p = len(pairs)
	for f in FACTORS:
		cg = CoarseGrid(native, native.res * f)
		out['resolutions'].append({'res': cg.res, 'factor': f, 'width': cg.width, 'height': cg.height,
		                           'memory_mb': memory_mb(cg, p)})
	cg = CoarseGrid(native, s['res'] or native.res)
	cg_ref = CoarseGrid(inp['ref_native'], cg.res) if inp['ref'] else cg
	out.update({'res': cg.res, 'factor': cg.factor, 'grid': {'width': cg.width, 'height': cg.height},
	            'memory_mb': memory_mb(cg, p)})
	eco = None
	try:
		eco = ecoregion_features(inp['eco_file']) if inp['eco_file'] else None
		out['ecoregions'] = {k: (rel(v) if k == 'file' else v) for k, v in eco.items()} if eco else None
		if not eco:
			out['warnings'].append('no file of the ecoregions: the ecoregion variables are taken as one ecoregion '
			                       '(TerrEco = 1)')
		elif inp['area']['kind'] == 'ecoregions' and len(eco['features']) > 1:
			_layer_in(eco['file'], cg.srs(), cg.extent)  # the check of run(): can they be drawn on the grid?
	except (ValueError, RuntimeError) as e:
		out['errors'].append(f'the file of the ecoregions cannot be used: {e}')
	local_res = local_resolution(pairs)
	for t in find_targets(inp['root'], s['cut'], s['ids']):
		try:
			vals, cov, hft_ids = hft_values(t, inp, cg, cg_ref, s['min_cover'], local_res)
		except ValueError as e:
			out['errors'].append(str(e))
			continue
		areas = hft_areas(t['shp'])
		hfts = []
		for h in hft_ids:
			n = int(len(vals[h]))
			st = 'too_few_pixels' if n < out['min_pixels'] else 'few_pixels' if n < out['few_pixels'] else 'ok'
			hfts.append({'hft': h, 'area_ha': round(areas.get(h, 0.0), 1), 'ref_pixels': n, 'status': st})
		folder = _run_folder(t, cg, inp['area']['tag'])
		other = _other_area(folder, inp)
		if other:
			out['errors'].append(f"{rel(folder)} has the results of another area ({other}), whose folder has the same "
			                     f"name: rename the folder of this area")
		out['targets'].append({'id': t['id'], 'cut': t['cut'], 'suggested': t['suggested'], 'shp': rel(t['shp']),
		                       'folder': rel(folder), 'exists': os.path.isfile(os.path.join(folder, RUN_JSON)),
		                       'in_area': cov is not None, 'hfts': hfts})
		if hfts and not any(x['ref_pixels'] for x in hfts):
			out['errors'].append(
				f"ID {t['id']}: the HFTs have no pixels with data in {rel(inp['ref'] or inp['area']['folder'])}: the "
				f"area does not contain the study area, and there is no {STUDY_AT_AREA} folder in it with the study area "
				f"at its pixel (get_vars_gee.py downloads it with the area)")
			continue
		skipped = [str(x['hft']) for x in hfts if x['status'] == 'too_few_pixels']
		if skipped:
			out['warnings'].append(f"ID {t['id']}: HFTs with too few pixels of {fmt_res(cg.res)} m, which will not be "
			                       f"compared: {', '.join(skipped)}")
	return out


def run(s, cfg, plots=True):
	t_start = time.time()
	gdal.SetCacheMax(512 * 2 ** 20)  # the block cache of GDAL (5 % of the memory by default) adds to the arrays
	inp = _inputs(s, cfg)
	root, pairs, source, area = inp['root'], inp['pairs'], inp['source'], inp['area']
	if inp['missing']:
		raise ValueError(f"these variables of the HFTs have no file in {rel(area['folder'])} (same name without the "
		                 f"number in front): " + ', '.join(inp['missing']))
	native, nodata = inp['native'], inp['nodata']
	cg = CoarseGrid(native, s['res'] or native.res)
	cg_ref = CoarseGrid(inp['ref_native'], cg.res) if inp['ref'] else cg
	targets = find_targets(root, s['cut'], s['ids'])
	eco = ecoregion_features(inp['eco_file']) if inp['eco_file'] else None
	marine = user_path(s['marine']) if str(s['marine'] or '').strip() else None
	if marine and not os.path.isfile(marine):
		raise ValueError(f'the file of the marine ecoregions does not exist: {marine}')
	p = len(pairs)
	print(f"Study area {s['study_area']}, scenario {s['scenario']} ({len(targets)} ID(s)); variables from {source}:")
	for pr in pairs:
		print(f"  {pr['key']}: {rel(pr['local'])} <-> {rel(pr['ecoregion'])}")
	what = 'Comparison area' if area['kind'] == 'comparison area' else 'Ecoregions'
	print(f"{what}: {rel(area['folder'])}{' (' + area['name'] + ')' if area['name'] else ''}")
	if inp['ref']:
		print(f"Statistics of the HFTs from {rel(inp['ref'])} (the study area at the pixel of the area)")
	elif inp['ref_note']:
		print(f"WARNING: {rel(os.path.join(area['folder'], STUDY_AT_AREA))} is not used: {inp['ref_note']}")
	print(f'Grid of the similarity: {fmt_res(cg.res)} m ({cg.factor} x {cg.factor} pixels of {fmt_res(native.res)} m), '
	      f'{cg.width} x {cg.height} pixels, {cg.crs_label()}')
	if eco:
		ids = [str(f['id']) for f in eco['features']]
		print(f"Ecoregions of the THDI: {rel(eco['file'])} ({', '.join(ids[:10])}{' ...' if len(ids) > 10 else ''}"
		      f"{f'; {len(ids)} polygons' if len(ids) > 10 else ''})")
	else:
		print('WARNING: no file of the ecoregions: the ecoregion variables are taken as one ecoregion (TerrEco = 1)')
	print('Reading the variables of the area...', flush=True)
	t0 = time.time()
	valid, X = load_stack([pr['ecoregion'] for pr in pairs], nodata, cg)
	if not len(X):
		raise ValueError(f"the variables of {rel(area['folder'])} have no pixel with data in all of them")
	eco_sd = np.array([X[:, j].std(dtype=np.float64) for j in range(p)])  # column by column: less memory
	print(f'  {len(X):,} pixels with data in the {p} variables ({time.time() - t0:.0f} s)', flush=True)
	flat = [pairs[j]['key'] for j in range(p) if not eco_sd[j] > 0]
	if flat:
		print(f"WARNING: constant in the whole area, so left out of the distance of every HFT: {', '.join(flat)}",
		      flush=True)
	# the ecoregion of every HFT, only when the area holds several ecoregions
	restrict = area['kind'] == 'ecoregions' and eco is not None and len(eco['features']) > 1
	eco_v = ecoregion_ids(eco, cg, valid) if restrict else None
	single = eco['features'][0]['id'] if area['kind'] == 'ecoregions' and eco and len(eco['features']) == 1 else None
	from scipy import ndimage
	ctx = {'X': X, 'valid': valid, 'cg': cg, 'pairs': pairs, 'eco_sd': eco_sd, 'stat': s['stat'],
	       'sim_label': 'MedianHFTSim' if s['stat'] == 'median' else 'AveHFTSim',
	       'structure': ndimage.generate_binary_structure(2, 2 if s['neighbours'] == 8 else 1),
	       'min_pixels': max(5, p + 1), 'few_pixels': 3 * p}
	local_res = local_resolution(pairs)
	by_dir = {}
	for t in targets:
		by_dir.setdefault(t['dir'], []).append(t)
	for d, ts in by_dir.items():  # before anything is replaced: the area of the folder and the HFTs of every ID
		other = _other_area(_run_folder(ts[0], cg, area['tag']), inp)
		if other:
			raise ValueError(f"{rel(_run_folder(ts[0], cg, area['tag']))} has the results of another area ({other}), "
			                 f"whose folder has the same name: rename the folder of this area")
		for t in ts:
			hft_values(t, inp, cg, cg_ref, s['min_cover'], local_res, valid, X)
	folders = []
	for d, ts in by_dir.items():
		outdir = _run_folder(ts[0], cg, area['tag'])
		old = _read_json(os.path.join(outdir, RUN_JSON))
		if s['ids'] and old.get('ids'):
			# the folder of a cut, pixel size and area keeps all its IDs: the others are computed again too
			others = [i for i in old['ids'] if i not in {t['id'] for t in ts}]
			again = []
			for i in others:
				try:
					again += [x for x in find_targets(root, ts[0]['cut'], [i]) if x['dir'] == d]
				except ValueError:
					pass
			if again:
				print(f"{rel(outdir)} also has the results of the ID(s) {', '.join(x['id'] for x in again)}: they are "
				      f"computed again too, so that the folder keeps all its IDs", flush=True)
				ts = ts + again
		final, outdir = outdir, outdir + '.tmp'  # written aside: a run that fails keeps the previous results
		if os.path.isdir(outdir):
			shutil.rmtree(outdir)
		os.makedirs(outdir)
		rows, recs = [], []
		for t in ts:
			vals, cov, hft_ids = hft_values(t, inp, cg, cg_ref, s['min_cover'], local_res, valid, X)
			ctx['areas'] = hft_areas(t['shp'])
			print(f"ID {t['id']}, cut {t['cut']}{' (suggested)' if t['suggested'] else ''}: {len(hft_ids)} HFTs "
			      f"({os.path.basename(t['shp'])})", flush=True)
			if cov is None:
				print('  (the HFTs are outside the area: all its pixels are compared)', flush=True)
			for h in hft_ids:
				cells = _cells(cov, h, cg)
				eco_of_hft = single
				if eco_v is not None:
					ok = cells[valid.ravel()[cells]] if len(cells) else cells
					rr, cc = ok // cg.width, ok % cg.width
					ids_h = eco_v[_compact_index(valid, rr, cc)] if len(ok) else np.zeros(0, dtype=np.int32)
					ids_h = ids_h[ids_h != NO_ECO]
					if ids_h.size:
						u, c = np.unique(ids_h, return_counts=True)
						eco_of_hft = int(u[c.argmax()])
					else:
						eco_of_hft = None
				row, rec = run_hft(h, ctx, t, vals[h], cells, eco_v, eco_of_hft, outdir)
				rows.append(row)
				recs.append(rec)
		hri = os.path.join(outdir, HRI_CSV)
		write_hri(hri, rows, pairs, ctx['sim_label'])
		thdi = os.path.join(outdir, THDI_CSV)
		cmd = ['Rscript', THDI_R, f'--hri={hri}', f'--out={thdi}'] + [f"--hfts={t['id']}={t['shp']}" for t in ts]
		if eco:
			cmd += [f"--ecoregions={eco['file']}", f"--eco-field={eco['field']}"]
		if marine:
			cmd += [f'--marine={marine}']
		print('Indices of the study area (thdi.R)...', flush=True)
		_r(cmd)
		plot_files = []
		if plots:
			cmd = ['Rscript', PLOTS_R, f'--dir={outdir}'] + [f"--study={t['id']}={t['shp']}" for t in ts]
			cmd += [f"--colours={t['id']}={t['colours']}" for t in ts if t['colours']]
			outline = area['polygons'] if area['tag'] else (eco['file'] if eco else None)
			if outline:
				cmd += [f'--ecoregions={outline}']
			if area['tag']:
				cmd += ['--area=' + (f"the comparison area {area['name']}" if area['kind'] == 'comparison area' and
				                     area['name'] else 'the comparison area' if area['kind'] == 'comparison area'
				                     else f"{area['tag']}")]
			print('Plots (similarity_plots.R)...', flush=True)
			if _r(cmd, required=False):
				plot_files = sorted(f for f in os.listdir(outdir) if f.lower().endswith('.png'))
		record = {'version': 1, 'date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'study_area': s['study_area'],
		          'scenario': s['scenario'], 'cut': ts[0]['cut'], 'ids': [t['id'] for t in ts], 'res': cg.res,
		          'native_res': native.res, 'factor': cg.factor, 'crs': cg.crs_label(),
		          'grid': {'width': cg.width, 'height': cg.height, 'geotransform': list(cg.gt)},
		          'valid_cells': int(len(X)), 'area': _area_record(inp),
		          'options': {'stat': s['stat'], 'neighbours': s['neighbours'], 'min_cover': s['min_cover'],
		                      'p_min': P_MIN, 'min_valid_fraction': MIN_VALID_FRACTION},
		          'variables_from': source,
		          'variables': [{'name': pr['name'], 'key': pr['key'], 'local': rel(pr['local']),
		                         'ecoregion': rel(pr['ecoregion'])} for pr in pairs],
		          'ecoregions': {'file': rel(eco['file']), 'field': eco['field'], 'features': eco['features']} if eco
		          else None,
		          'marine': rel(marine),
		          'shapefiles': {t['id']: rel(t['shp']) for t in ts},
		          'colours': {t['id']: rel(t['colours']) for t in ts},
		          'hfts': recs, 'tables': {'hri': HRI_CSV, 'thdi': THDI_CSV if os.path.isfile(thdi) else None},
		          'plots': plot_files, 'seconds': round(time.time() - t_start, 1)}
		write_metadata(os.path.join(outdir, RUN_TXT), record, thdi)
		with open(os.path.join(outdir, RUN_JSON), 'w', encoding='utf-8') as f:
			json.dump(record, f, indent=1, ensure_ascii=False)
		if os.path.isdir(final):
			shutil.rmtree(final)  # a previous run with the same cut, pixel size and area: replaced
		os.replace(outdir, final)
		folders.append(final)
	return folders


def _read_json(path):
	try:
		with open(path, encoding='utf-8') as f:
			return json.load(f)
	except (OSError, ValueError):
		return {}


def _r(cmd, required=True):
	"""Runs an R script; stops (or, if not REQUIRED, warns) when it fails."""
	print('>>> ' + ' '.join(os.path.basename(c) if i < 2 else c for i, c in enumerate(cmd)), flush=True)
	try:
		r = subprocess.run(cmd, cwd=HERE)
	except OSError as e:
		r = None
		err = f'{cmd[1]} could not be run ({e}): is R (Rscript) installed in the environment?'
	if r is not None and r.returncode == 0:
		return True
	err = err if r is None else f'{os.path.basename(cmd[1])} failed (exit code {r.returncode}): see the messages above'
	if required:
		sys.exit('ERROR: ' + err)
	print('WARNING: ' + err, flush=True)
	return False


def write_metadata(path, rec, thdi):
	lines = ['eHabitat+ similarity metadata: options and inputs of the files in this folder (similarity.py)',
	         f"Date: {rec['date']}", f"Study area: {rec['study_area']}", f"Scenario: {rec['scenario']}",
	         f"Cut: {rec['cut']}", f"IDs: {', '.join(rec['ids'])}", '',
	         '## Options',
	         f"ECOSIM_RESOLUTION = {fmt_res(rec['res'])} m ({rec['factor']} x {rec['factor']} pixels of "
	         f"{fmt_res(rec['native_res'])} m)",
	         f"ECOSIM_STAT = {rec['options']['stat']!r}", f"ECOSIM_NEIGHBOURS = {rec['options']['neighbours']}",
	         f"ECOSIM_MIN_COVER = {rec['options']['min_cover']}",
	         f"Similarity set to 0 at or below p = {P_MIN}; a coarser pixel needs data in "
	         f"{int(MIN_VALID_FRACTION * 100)} % of its area", '',
	         '## Grid of the similarity',
	         f"{rec['crs']}, {rec['grid']['width']} x {rec['grid']['height']} pixels, geotransform "
	         f"{rec['grid']['geotransform']}", f"Pixels with data in all the variables: {rec['valid_cells']}", '',
	         f"## Variables (from {rec['variables_from']}): local <-> ecoregion"]
	lines += [f"{v['key']}: {v['local']} <-> {v['ecoregion']}" for v in rec['variables']]
	a = rec['area']
	lines += ['', '## Area compared', f"{a['kind']}: {a['folder']}{' (' + a['name'] + ')' if a['name'] else ''}",
	          f"Statistics of the HFTs from: {a['reference']}", '', '## Ecoregions (ecoregion of each HFT and THDI)']
	if rec['ecoregions']:
		lines.append(f"{rec['ecoregions']['file']} (field {rec['ecoregions']['field']}): " + '; '.join(
			f"{f['id']} {f['name'] or ''}".strip() for f in rec['ecoregions']['features']))
	else:
		lines.append('none: the ecoregion variables are taken as one ecoregion')
	lines.append(f"Marine ecoregions: {rec['marine'] or 'none'}")
	lines += ['', '## HFTs (shapefiles: ' + ', '.join(rec['shapefiles'].values()) + ')',
	          'id hft pixels variables ridge threshold status left_out(constant)']
	for h in rec['hfts']:
		lines.append(f"{h['id']} {h['hft']} {h['ref_pixels']} {h['vars_used'] if h['vars_used'] is not None else 'NA'} "
		             f"{h['ridge'] if h['ridge'] is not None else 'NA'} "
		             f"{_num(h['threshold'])} {h['status']} {','.join(h.get('vars_dropped') or []) or '-'}")
	if os.path.isfile(thdi):
		lines += ['', '## Indices of the study areas (thdi_res.csv)'] + open(thdi).read().strip().splitlines()
	lines += ['', f"Time: {rec['seconds']} s"]
	with open(path, 'w', encoding='utf-8') as f:
		f.write('\n'.join(lines) + '\n')


def main():
	ap = argparse.ArgumentParser(description='Similarity of the HFTs in their ecoregion, landscape metrics and THDI.')
	ap.add_argument('--study-area', help='study area of the results (default: STUDY_AREA of config.py)')
	ap.add_argument('--scenario', help='scenario of the results (default: that of config.py)')
	ap.add_argument('--area', help="folder (in inVars) of the ecoregion or comparison area, '' = the folder ecoregion "
	                               "of the variables of the HFTs (default: ECOSIM_AREA)")
	ap.add_argument('--cut', help="cut of the tree: p10 ... p90, '' = the suggested one (default: ECOSIM_CUT)")
	ap.add_argument('--id', nargs='+', dest='ids', help='IDs of the study areas (default: all those with HFTs; the '
	                                                    'other IDs of a folder of results that already has them are '
	                                                    'computed again too)')
	ap.add_argument('--res', type=float, help='pixel size of the similarity in metres, 0 = that of the ecoregion '
	                                          'variables (default: ECOSIM_RESOLUTION)')
	ap.add_argument('--stat', choices=('mean', 'median'), help='threshold of the similar areas (default: ECOSIM_STAT)')
	ap.add_argument('--neighbours', type=int, choices=(4, 8), help='neighbours of the patches (ECOSIM_NEIGHBOURS)')
	ap.add_argument('--min-cover', type=float, help='minimum cover of a pixel by its HFT (ECOSIM_MIN_COVER)')
	ap.add_argument('--ecoregions', help='polygons of the ecoregions (ECOSIM_ECOREGIONS)')
	ap.add_argument('--marine', help='polygons of the marine ecoregions for the THDI (ECOSIM_MARINE)')
	ap.add_argument('--no-plots', action='store_true', help='do not run similarity_plots.R')
	ap.add_argument('--plan', action='store_true', help='only print what would be done (JSON)')
	a = ap.parse_args()
	cfg = hft_state.read_config()
	try:
		if a.plan:
			print(json.dumps(plan(a.study_area, a.scenario, a.cut, a.ids, a.res, a.min_cover, a.ecoregions, cfg,
			                      area=a.area), indent=1, ensure_ascii=False))
			return
		s = _settings(cfg, a.study_area, a.scenario, a.cut, a.ids, a.res, a.stat, a.neighbours, a.min_cover,
		              a.ecoregions, a.marine, a.area)
		folders = run(s, cfg, plots=not a.no_plots)
	except (ValueError, RuntimeError) as e:  # RuntimeError: GDAL/OGR (e.g. a file that cannot be read)
		sys.exit(f'ERROR: {e}')
	for f in folders:
		print(f'Results in {rel(f)}')


if __name__ == '__main__':
	os.chdir(HERE)
	main()
