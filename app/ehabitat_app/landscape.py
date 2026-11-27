"""Landscape metrics of maps of patches of the user (segm_optim/landscape_metrics.py, without the similarity): the
catalogue of the metrics (segm_optim/lsm.py, those of the R package landscapemetrics), what a raster or vector file
holds, the runs (a job, in segm_optim/landscape/<name>) and their results: tables, plots and the patches of each map,
as GeoJSON for the map of the app."""

import heapq
import importlib
import json
import os
import re
import shutil
import sys
import threading

from . import inputs, jobs, paths, results

ENGINE = os.path.join(paths.SEGM, 'landscape_metrics.py')
OUT = os.path.join(paths.SEGM, 'landscape')
RUN_JSON = 'landscape_run.json'
NAME_RE = re.compile(r'^[A-Za-z0-9_-]{1,80}$')
MAX_PATCHES = 40000  # polygons of a map drawn on the map of the app: the largest ones beyond them
PLOTS = [(r'^metrics_trend\.png$', 'Class metrics in every map'),
         (r'^landscape_trend\.png$', 'Landscape-level metrics in every map')]
_lock = threading.Lock()
_mtimes = {}


def _module(name):
	"""The module NAME of segm_optim (lsm or landscape_metrics), imported again when its file changed (an update)."""
	with _lock:
		path = os.path.join(paths.SEGM, name + '.py')
		if not os.path.isfile(path):
			raise ValueError(f'this copy of eHabpy has no segm_optim/{name}.py: update it (Check for updates)')
		mtime = os.path.getmtime(path)
		try:
			mod = importlib.import_module(name)
			if _mtimes.get(name) not in (None, mtime):
				if name == 'landscape_metrics':
					importlib.reload(importlib.import_module('lsm'))
				mod = importlib.reload(mod)
		except (ImportError, SyntaxError) as e:
			raise ValueError(f'segm_optim/{name}.py cannot be loaded: {e}')
		_mtimes[name] = mtime
		return mod


def _plain(x):
	if isinstance(x, dict):
		return {str(k): _plain(v) for k, v in x.items()}
	if isinstance(x, (list, tuple)):
		return [_plain(v) for v in x]
	if hasattr(x, 'item') and callable(x.item) and not isinstance(x, (str, bytes)):
		try:
			return x.item()
		except (TypeError, ValueError):
			return str(x)
	return x


def catalogue():
	"""The metrics (key, level, name, type, unit, what they measure, how to read them, whether larger values mean more
	fragmentation, their page in the manual of landscapemetrics), the default set of a run and the references; and the
	columns of hri_results.csv of the similarity with the metrics they are."""
	lsm = _module('lsm')
	eng = _module('landscape_metrics')
	sim = []
	try:
		import similarity
		sim = [{'column': c, 'key': 'c_' + m} for c, m in getattr(similarity, 'SIM_LSM', [])]
	except Exception:  # noqa: BLE001 - the similarity is optional here
		pass
	return {'metrics': [dict(m, docs=lsm.DOCS.format(function=m['function'])) for m in lsm.CATALOGUE],
	        'fragmentation': list(eng.FRAGMENTATION), 'references': [lsm.REFERENCE, lsm.FRAGSTATS], 'similarity': sim}


def inspect(path, layer=None):
	return _plain(_module('landscape_metrics').inspect(path, layer))


def _folder(name):
	name = str(name or '')
	if not NAME_RE.match(name):
		raise ValueError(f'{name!r} is not the name of a run')
	d = os.path.join(OUT, name)
	if not paths.inside(d, OUT):
		raise ValueError(f'{name!r} is not the name of a run')
	return d


def runs():
	out = []
	for r in _module('landscape_metrics').runs():
		out.append({k: r[k] for k in ('name', 'date', 'mode', 'maps')} | {'folder': paths.rel(r['folder'])})
	return out


def run(spec):
	"""Run landscape_metrics.py with SPEC (a job); its result is the name of the run."""
	if jobs.running('landscape'):
		raise ValueError('the landscape metrics are already being computed')
	eng = _module('landscape_metrics')
	spec = dict(spec or {})
	name = eng._check_spec(spec)[0]
	_folder(name)
	os.makedirs(OUT, exist_ok=True)
	tmp = os.path.join(OUT, f'.spec_{name}.json')
	with open(tmp, 'w', encoding='utf-8') as f:
		json.dump(spec, f, indent=1, ensure_ascii=False)
	cmd = [sys.executable, 'landscape_metrics.py', '--spec', tmp]

	def work(job):
		try:
			code = jobs.run_command(job, cmd, cwd=paths.SEGM)
			if code != 0:
				raise RuntimeError(f'landscape_metrics.py failed (exit code {code}): see the messages above')
			return {'name': name}
		finally:  # also when it fails or is stopped: its spec and the folder it was writing
			try:
				os.remove(tmp)
			except OSError:
				pass
			shutil.rmtree(os.path.join(OUT, name + '.tmp'), ignore_errors=True)

	job = jobs.start_function('landscape', 'Landscape metrics', work)
	return job.id


def _rows(table):
	if not table:
		return []
	h = table['header']
	return [{k: (v if v != 'NA' else None) for k, v in zip(h, r)} for r in table['rows']]


def run_results(name):
	"""A finished run: its record, the tables of the classes and the landscape (all their rows) and of the patches
	(the first rows), the plots and landscape_metadata.txt."""
	d = _folder(name)
	p = os.path.join(d, RUN_JSON)
	if not os.path.isfile(p):
		raise ValueError(f'{paths.rel(d)} has no finished run ({RUN_JSON})')
	with open(p, encoding='utf-8') as f:
		rec = json.load(f)
	tables = {}
	for key, limit in (('class', 100000), ('landscape', 100000), ('patch', 2000)):
		f = (rec.get('tables') or {}).get(key)
		t = os.path.join(d, f) if f else None
		tables[key] = results._table(t, limit=limit) if t and os.path.isfile(t) else None
		if tables[key]:
			with open(t, encoding='utf-8', errors='replace') as fh:
				tables[key]['total'] = max(0, sum(1 for _ in fh) - 1)
	plots = []
	for f in rec.get('plots') or []:
		fp = os.path.join(d, f)
		if os.path.isfile(fp):
			label = next((lab for pat, lab in PLOTS if re.match(pat, f)), f)
			plots.append({'label': label, 'file': f, 'url': results._url(fp)})
	meta = os.path.join(d, 'landscape_metadata.txt')
	return _plain({'name': name, 'folder': paths.rel(d), 'path': d, 'run': rec, 'tables': tables,
	               'class': _rows(tables['class']), 'landscape': _rows(tables['landscape']), 'plots': plots,
	               'metadata': open(meta, encoding='utf-8', errors='replace').read() if os.path.isfile(meta) else None})


def patches(name, index):
	"""The patches of the map INDEX (1, 2...) of a run as GeoJSON in WGS 84 (simplified to their pixels), the largest
	MAX_PATCHES of them, with their class and patch metrics; and the extent of the map."""
	from osgeo import ogr, osr
	d = _folder(name)
	with open(os.path.join(d, RUN_JSON), encoding='utf-8') as f:
		rec = json.load(f)
	maps = rec.get('maps') or []
	i = int(index)
	if not 1 <= i <= len(maps):
		raise ValueError(f'the run has no map {index}')
	m = maps[i - 1]
	gp = os.path.join(d, os.path.basename(m.get('patches_gpkg') or f'patches_{i}.gpkg'))
	if not os.path.isfile(gp):
		raise ValueError(f'{paths.rel(gp)} not found')
	ds = ogr.Open(gp)
	lyr = ds.GetLayer(0)
	srs = lyr.GetSpatialRef().Clone()
	srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
	ct = osr.CoordinateTransformation(srs, inputs._wgs84())
	defn = lyr.GetLayerDefn()
	names = [defn.GetFieldDefn(k).GetName() for k in range(defn.GetFieldCount())]
	heap, total, seq = [], 0, 0
	for feat in lyr:
		g = feat.GetGeometryRef()
		if g is None:
			continue
		total += 1
		seq += 1
		area = feat.GetField('area') if 'area' in names else 0
		item = (area or 0, seq, feat.GetFID())
		if len(heap) < MAX_PATCHES:
			heapq.heappush(heap, item)
		elif item[0] > heap[0][0]:
			heapq.heapreplace(heap, item)
	res = (m.get('grid') or {}).get('res') or 0
	feats = []
	for _, _, fid in sorted(heap, key=lambda x: x[1]):
		feat = lyr.GetFeature(fid)
		g = feat.GetGeometryRef().Clone()
		if res:
			g = g.SimplifyPreserveTopology(res / 2)
		g.Transform(ct)
		props = {k: feat.GetField(k) for k in names}
		feats.append({'type': 'Feature', 'properties': props, 'geometry': inputs._json(g)})
	xmin, xmax, ymin, ymax = lyr.GetExtent()
	ext = (m.get('grid') or {}).get('extent')
	bbox = inputs._bbox_wgs84(srs, *(ext if ext else (xmin, ymin, xmax, ymax)))
	return _plain({'map': m.get('map'), 'index': i, 'patches': {'type': 'FeatureCollection', 'features': feats},
	               'total': total, 'truncated': total > MAX_PATCHES, 'max_patches': MAX_PATCHES, 'bbox': bbox,
	               'classes': m.get('classes') or []})


def boundary_geojson(name):
	"""The boundary of a run (GeoJSON in WGS 84), or None."""
	d = _folder(name)
	with open(os.path.join(d, RUN_JSON), encoding='utf-8') as f:
		rec = json.load(f)
	b = rec.get('boundary')
	if not b or not b.get('file'):
		return None
	p = b['file'] if os.path.isabs(b['file']) else os.path.join(paths.ROOT, b['file'])
	if not os.path.isfile(p):
		return None
	return inputs.vector_geojson(p, b.get('layer'))


def remove(name):
	"""Delete the folder of a run."""
	d = _folder(name)
	if jobs.running('landscape'):
		raise ValueError('wait until the landscape metrics being computed finish')
	if not os.path.isdir(d):
		raise ValueError(f'no run {name}')
	shutil.rmtree(d)
	return True
