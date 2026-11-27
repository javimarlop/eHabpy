#### eHabitat+: test of the landscape metrics (lsm.py) against the R package landscapemetrics
#### License: CC BY-SA 3.0
"""
Checks that lsm.py gives the values of landscapemetrics 2.2.1 (calculate_lsm) for every metric of its catalogue, at
the patch, class and landscape levels, in three small maps (four classes with cells without data, a habitat map of 0
and 1, and scattered classes split by a row without data), with patches of 8 neighbours and an edge depth of 1 cell,
and of 4 neighbours and 2 cells (lsm_landscapemetrics.json, made with landscapemetrics in R); and a few values worked
out by hand.

The only metric that may differ is CIRCLE: lsm.py takes the smallest circle that contains every corner of the cells of
a patch, while landscapemetrics 2.2.1 gives a slightly smaller circle for a few patches (one that does not contain the
whole patch); the maps of this test have none of them.

    cd segm_optim && python -s tests/test_lsm.py
"""

import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import lsm  # noqa: E402


def same(a, b):
	if a is None or b is None:
		return a is None and b is None
	return abs(a - b) <= 1e-7 * max(1.0, abs(a), abs(b))


def test_landscapemetrics():
	ref = json.load(open(os.path.join(HERE, 'lsm_landscapemetrics.json'), encoding='utf-8'))
	checked, bad = 0, []
	for key, exp in ref['runs'].items():
		name, d, e = key.rsplit('_', 2)
		m = np.array(ref['maps'][name]['map'])
		out = lsm.compute(m, m != -1, ref['maps'][name]['res'], directions=int(d), edge_depth=int(e), classes_max=6)
		classes = {str(r['class']): r for r in out['class']}
		for metric, by_class in exp['class'].items():
			for c, v in by_class.items():
				mine = classes[c][metric] if c in classes else (None if all(x[metric] is None for x in classes.values())
				                                                 else 'not NA')
				checked += 1
				if mine == 'not NA' or not same(mine, v):
					bad.append(f'{key} class {c} {metric}: {mine} != {v}')
		for metric, v in exp['landscape'].items():
			checked += 1
			if not same(out['landscape'][metric], v):
				bad.append(f"{key} landscape {metric}: {out['landscape'][metric]} != {v}")
		for k, vals in exp['patch'].items():
			c, metric = k.split(':')
			mine = sorted((p[metric] for p in out['patch'] if str(p['class']) == c), key=lambda z: (z is None, z or 0))
			checked += 1
			if len(mine) != len(vals) or not all(same(x, y) for x, y in zip(mine, vals)):
				bad.append(f'{key} patches of class {c} {metric}: {mine[:5]} != {vals[:5]}')
	assert not bad, f'{len(bad)} of {checked} values differ from landscapemetrics:\n' + '\n'.join(bad[:30])
	# every metric of the catalogue was compared
	assert {f'{lvl[0]}_{m}' for run in ref['runs'].values() for lvl in ('class', 'landscape') for m in run[lvl]} | \
		{f"p_{k.split(':')[1]}" for run in ref['runs'].values() for k in run['patch']} == set(lsm.keys())
	print(f'landscapemetrics: {checked} values, all equal')


def test_by_hand():
	# a 4 x 4 map: a 2 x 2 block of class 1 in a corner, the rest class 0; cells of 10 m
	m = np.zeros((4, 4), dtype=int)
	m[:2, :2] = 1
	out = lsm.compute(m, np.ones(m.shape, bool), 10)
	c1 = next(r for r in out['class'] if r['class'] == 1)
	assert same(c1['ca'], 4 * 100 / 1e4) and same(c1['pland'], 25) and c1['np'] == 1
	assert same(c1['te'], 4 * 10)  # its edges with class 0 (not those with the edge of the map)
	assert same(c1['ed'], 40 / 0.16) and same(c1['lsi'], 1) and same(c1['ai'], 100)
	assert same(c1['division'], 1 - 0.25 ** 2) and same(c1['mesh'], 0.04 ** 2 / 0.16) and same(c1['split'], 16)
	p1 = next(p for p in out['patch'] if p['class'] == 1)
	assert same(p1['perim'], 80) and same(p1['shape'], 1) and same(p1['core'], 0)  # all its cells are edge
	assert same(p1['gyrate'], math.sqrt(0.5) * 10) and same(p1['circle'], 1 - 400 / (math.pi * (20 * math.sqrt(2) / 2) ** 2))
	assert out['landscape']['pr'] == 2 and same(out['landscape']['shdi'], -(0.25 * math.log(.25) + .75 * math.log(.75)))
	# patches that only touch at a corner: one patch with 8 neighbours, two with 4
	m = np.array([[1, 0], [0, 1]])
	assert lsm.compute(m, np.ones(m.shape, bool), 1, ['c_np'], directions=8)['class'][1]['np'] == 1
	assert lsm.compute(m, np.ones(m.shape, bool), 1, ['c_np'], directions=4)['class'][1]['np'] == 2
	# ENN: two single cells of class 1 three cells apart (centre to centre); NA with one patch
	m = np.zeros((1, 7), dtype=int)
	m[0, [1, 4]] = 1
	out = lsm.compute(m, np.ones(m.shape, bool), 2, ['p_enn', 'c_enn_mn'])
	assert [p['enn'] for p in out['patch'] if p['class'] == 1] == [6.0, 6.0]
	assert same(out['class'][0]['enn_mn'], 4.0)  # class 0: columns 0, 2-3 and 5-6, two cells (4 m) apart
	m = np.array([[1, 1, 0]])
	assert lsm.compute(m, np.ones(m.shape, bool), 1, ['c_enn_mn'])['class'][1]['enn_mn'] is None  # one patch: NA
	print('by hand: ok')


def test_only_classes():
	"""The metrics of one class (as similarity.py computes those of the similar areas) are those of the whole map."""
	rng = np.random.default_rng(5)
	m = rng.integers(0, 3, (30, 40))
	valid = rng.random(m.shape) > 0.1
	keys = [k for k in lsm.keys('class')]
	a = lsm.compute(m, valid, 25, keys)
	b = lsm.compute(m, valid, 25, keys, only_classes=[1])
	ra = next(r for r in a['class'] if r['class'] == 1)
	assert len(b['class']) == 1 and not b['landscape']
	for k, v in b['class'][0].items():
		assert same(v, ra[k]) if v is not None else ra[k] is None, k
	print('one class: ok')


if __name__ == '__main__':
	test_landscapemetrics()
	test_by_hand()
	test_only_classes()
	print('OK')
