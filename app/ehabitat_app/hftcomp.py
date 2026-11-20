"""HFT comparator: the tables of the HFTs (hft_means_...csv, written by moranvar_plots.R in the folder of every scenario and
cut) of several runs, put together in the folder segm_optim/comparison, and the run of hft_compare.R, which joins them and
draws the NMDS, the tree, the PCA and the similarity between all the HFTs (results in comparison/results)."""

import csv
import os
import re
import shutil

from . import jobs, paths, results

DIR = os.path.join(paths.SEGM, 'comparison')
OUT = os.path.join(DIR, 'results')
IDS = ('study_area', 'scenario', 'id', 'threshold', 'cut', 'suggested')
PLOTS = [
	(r'^nmds_hfts\.png$', 'nmds', 'NMDS of the HFTs'),
	(r'^tree_hfts\.png$', 'tree', 'Tree of the HFTs (Ward)'),
	(r'^heatmap_hfts\.png$', 'heatmap', 'Profiles of the HFTs (standardized means)'),
	(r'^pca_biplot_hfts\.png$', 'pca', 'PCA of the HFTs'),
	(r'^similarity_hfts\.png$', 'similarity', 'Similarity between the HFTs'),
	(r'^similarity_sources\.png$', 'sources', 'Similarity between the sources'),
	(r'^radarplot_hfts_(?P<s>.+)\.png$', 'radar', 'Radar plot of the source {s}'),
]
TABLES = [('hfts_all.csv', 'The HFTs of all the tables, with their unique id (label) and code'),
          ('best_matches.csv', 'The most similar HFT of every other source'),
          ('similarity_sources.csv', 'Mean similarity between sources'),
          ('nmds_hfts_variables.csv', 'Variables in the NMDS (envfit)'),
          ('pca_biplot_hfts_variance.csv', 'Variance of each axis of the PCA')]


def _safe(text):
	return re.sub(r'[^A-Za-z0-9_-]', '_', str(text))


def _read(path, limit=None):
	"""(header, rows) of a table separated by spaces, commas or semicolons."""
	with open(path, newline='', encoding='utf-8-sig', errors='replace') as f:
		first = f.readline()
		f.seek(0)
		delim = ';' if ';' in first else ',' if ',' in first else ' '
		rows = []
		for r in csv.reader(f, delimiter=delim):
			if r:
				rows.append(r)
			if limit and len(rows) > limit:
				break
	return (rows[0], rows[1:]) if rows else ([], [])


def _info(path):
	"""What a CSV of the comparison folder is: a table of HFTs (label and area_pct) with its sources, or not."""
	rec = {'file': os.path.basename(path), 'size': os.path.getsize(path), 'mtime': os.path.getmtime(path),
	       'usable': False, 'reason': None, 'hfts': 0, 'sources': [], 'variables': 0, 'resolution': None}
	try:
		header, rows = _read(path)
	except OSError as e:
		rec['reason'] = str(e)
		return rec
	if 'label' not in header or 'area_pct' not in header:
		rec['reason'] = 'not a table of HFTs (hft_means_...): it has no label and area_pct columns'
		return rec
	ix = {h: i for i, h in enumerate(header)}
	rec.update(usable=True, hfts=len(rows), variables=len(header) - ix['area_pct'] - 1)
	src = {}
	for r in rows:
		label = r[ix['label']] if ix['label'] < len(r) else ''
		parts = label.split('/')
		key = '/'.join(parts[:-1]) if len(parts) > 1 else label
		src[key] = src.get(key, 0) + 1
		if 'resolution' in ix and ix['resolution'] < len(r) and rec['resolution'] is None:
			rec['resolution'] = r[ix['resolution']]
	rec['sources'] = [{'source': k, 'hfts': v} for k, v in src.items()]
	return rec


def files():
	"""The CSV files of the comparison folder, with what each one has."""
	os.makedirs(DIR, exist_ok=True)
	out = [_info(os.path.join(DIR, f)) for f in sorted(os.listdir(DIR)) if f.lower().endswith('.csv')
	       and os.path.isfile(os.path.join(DIR, f))]
	res = {r['resolution'] for r in out if r['usable'] and r['resolution'] not in (None, 'NA')}
	return {'folder': paths.rel(DIR), 'path': DIR, 'files': out, 'resolutions': sorted(res),
	        'results': os.path.isfile(os.path.join(OUT, 'hfts_all.csv'))}


def _scenario_dir(study_area, scenario):
	d = os.path.join(paths.RESULTS, study_area, scenario)
	if not paths.inside(d, paths.RESULTS) or not os.path.isdir(d):
		raise ValueError(f'{study_area}/{scenario} has no results')
	return d


def _means_in(folder):
	if not os.path.isdir(folder):
		return []
	return sorted(f for f in os.listdir(folder) if re.match(r'^hft_means_.*\.csv$', f))


def sources(study_area, scenario):
	"""The tables of HFTs that a scenario has: the one of all the cuts and one per cut (with SIMIL = 3)."""
	d = _scenario_dir(study_area, scenario)
	out = []
	allc = _means_in(d)
	if allc:
		out.append({'which': 'all', 'label': 'All the cuts in one table', 'files': allc})
	cuts = []
	for sub in sorted(os.listdir(d)):
		p = os.path.join(d, sub)
		if re.match(r'^p\d+$', sub) and os.path.isdir(p) and _means_in(p):
			cuts.append((int(sub[1:]), sub, _means_in(p)))
	for _, sub, fs in sorted(cuts):
		opt = any('_optim' in f for f in fs)
		out.append({'which': sub, 'label': f'Cut {sub}' + (' (suggested)' if opt else ''), 'files': fs, 'suggested': opt})
	return out


def _target(path):
	"""Name to copy a table to the comparison folder: hft_means_<area>_<scenario>_<ID>_<k>[_<cut>[_optim]].csv, built from
	its rows when the file has an older name (that of the run is not in it), so that tables of different runs never
	have the same name and need no renaming by hand."""
	name = os.path.basename(path)
	header, rows = _read(path, limit=2)
	if not rows or not set(IDS) <= set(header):
		return name
	ix = {h: i for i, h in enumerate(header)}
	first = rows[0]
	area, scen, ident = _safe(first[ix['study_area']]), first[ix['scenario']], first[ix['id']]
	if f'{area}_{scen}_' in name:
		return name
	kind = 'hft_stats' if name.startswith('hft_stats') else 'hft_means'
	try:
		k = str(int(round(float(first[ix['threshold']]) * 10)))
	except ValueError:
		k = '0'
	cuts = {r[ix['cut']] for r in _read(path)[1] if ix['cut'] < len(r)}
	tag = ''
	if len(cuts) == 1 and next(iter(cuts)) not in ('NA', ''):
		tag = '_' + next(iter(cuts)) + ('_optim' if first[ix['suggested']] == '1' else '')
	return f'{kind}_{area}_{scen}_{ident}_{k}{tag}.csv'


def add(study_area, scenario, which='suggested'):
	"""Copy the tables of HFTs of a scenario to the comparison folder. WHICH: 'suggested' (the table of the suggested cut,
	or the one of all the cuts when it has no cuts), 'all' (the table of all the cuts), 'each' (one per cut) or a cut (p50)."""
	d = _scenario_dir(study_area, scenario)
	avail = sources(study_area, scenario)
	if not avail:
		raise ValueError(f'{study_area}/{scenario} has no tables of HFTs (hft_means_...): run the model again (HFTs and plots)')
	if which == 'suggested':
		pick = [s for s in avail if s.get('suggested')] or [s for s in avail if s['which'] == 'all'] or avail[:1]
	elif which == 'each':
		pick = [s for s in avail if s['which'] != 'all'] or avail
	else:
		pick = [s for s in avail if s['which'] == which]
		if not pick:
			raise ValueError(f'{which} is not a table of {study_area}/{scenario}')
	os.makedirs(DIR, exist_ok=True)
	copied = []
	for s in pick:
		folder = d if s['which'] == 'all' else os.path.join(d, s['which'])
		for f in s['files']:
			src = os.path.join(folder, f)
			dst = os.path.join(DIR, _target(src))
			shutil.copyfile(src, dst)
			copied.append(os.path.basename(dst))
	return {'copied': copied}


def remove(name):
	p = os.path.join(DIR, os.path.basename(name))
	if not name.lower().endswith('.csv') or not os.path.isfile(p):
		raise ValueError(f'{name} is not a CSV file of the comparison folder')
	os.remove(p)
	return True


def run():
	if jobs.running('model'):
		raise ValueError('the model is running: compare the HFTs when it finishes')
	if jobs.running('hftcomp'):
		raise ValueError('the comparison is already running')
	usable = [f for f in files()['files'] if f['usable']]
	if not usable:
		raise ValueError('the comparison folder has no tables of HFTs (hft_means_...): add some first')
	if sum(f['hfts'] for f in usable) < 4:
		raise ValueError('at least 4 HFTs are needed to compare them')
	cmd = ['Rscript', 'hft_compare.R']

	def after(job):
		job.say(f'\nResults in {paths.rel(OUT)}'.replace(os.sep, '/'))
		return {'folder': paths.rel(OUT)}

	env = dict(os.environ, EHAB_PROGRESS='1')
	return jobs.start_commands('hftcomp', 'HFT comparator', [cmd], env=env, after=after).id


def run_results():
	"""The plots, tables and notes of the last comparison."""
	if not os.path.isfile(os.path.join(OUT, 'hfts_all.csv')):
		return {'path': OUT, 'plots': [], 'tables': [], 'metadata': None}
	names = sorted(os.listdir(OUT))
	plots = []
	for pat, kind, label in PLOTS:
		for f in names:
			m = re.match(pat, f)
			if m:
				plots.append({'kind': kind, 'label': label.format(**m.groupdict()), 'file': f, 'url': results._url(os.path.join(OUT, f))})
	tables = []
	for f, label in TABLES:
		p = os.path.join(OUT, f)
		if os.path.isfile(p):
			tables.append(dict(results._table(p, limit=3000), label=label))
	meta = os.path.join(OUT, 'comparison_metadata.txt')
	return {'path': OUT, 'folder': paths.rel(OUT), 'plots': plots, 'tables': tables,
	        'metadata': open(meta, encoding='utf-8', errors='replace').read() if os.path.isfile(meta) else None}
