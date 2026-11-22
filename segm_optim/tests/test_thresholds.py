#### eHabitat+ HFT model: test of the thresholds of the segmentation (SEGM_THRESHOLDS, hft_state.py)
#### License: CC BY-SA 3.0
"""
Checks the list of thresholds of the segmentation (SEGM_THRESHOLDS in config.py) as hft_state.py reads it: its
validation (sorted, without repeats), the index k of every threshold (its position in the list) and of a fixed
THRESHOLD_METHOD, the layout of segm_id (1 + 5 digits with the nine of 0.1-0.9, 2 + 7 with any other list) and its
largest value (segm_id_error), csv/segm_thresholds.csv and its fallback (the thresholds of config.py: k / 10 for the
segmentations of eHabpy 1.5 and before), the checks of THRESHOLD_METHOD and SMALL_SEGM_MERGE (options_error,
merge_rule), that the default values of SEGM_THRESHOLDS and SMALL_SEGM_MERGE are not recorded in the inputs of the
segmentation (the runs made before these options existed stay up to date), and that the app replaces a value of
config.py written over several lines whole (set_option).

    cd segm_optim && python -s tests/test_thresholds.py
"""

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim: hft_state.py

import hft_state as h  # noqa: E402

SMALL = [0.02, 0.04, 0.06, 0.08, 0.1]


def rejects(f, *args):
	try:
		f(*args)
	except ValueError as e:
		return str(e)
	raise AssertionError(f'{f.__name__}{args} did not raise ValueError')


def test_parse():
	assert h.parse_thresholds(SMALL) == [0.02, 0.04, 0.06, 0.08, 0.1]
	assert h.parse_thresholds(tuple(h.DEFAULT_THRESHOLDS)) == h.DEFAULT_THRESHOLDS
	assert len(h.parse_thresholds([i / 100 for i in range(1, 100)])) == h.MAX_THRESHOLDS
	for bad in ([], [0], [1], [0.5, 1.2], [-0.1], ['0.1'], [True], [None], [float('nan')], 0.5, '0.1, 0.2',
	            [i / 1000 for i in range(1, 101)]):
		assert rejects(h.parse_thresholds, bad).startswith('SEGM_THRESHOLDS: ')
	# not sorted, or repeated (also within THRESHOLD_TOL): an error, not changed silently
	assert 'from the lowest to the highest' in rejects(h.parse_thresholds, [0.4, 0.1])
	assert 'from the lowest to the highest' in rejects(h.parse_thresholds, list(reversed(h.DEFAULT_THRESHOLDS)))
	assert '0.1 is repeated' in rejects(h.parse_thresholds, [0.1, 0.1, 0.2])
	assert '0.3 is repeated' in rejects(h.parse_thresholds, [0.1, 0.3, 0.3 + 1e-12])
	assert h.segm_thresholds({}) == h.DEFAULT_THRESHOLDS  # a config.py without the option
	print('test_parse: OK')


def test_levels():
	assert h.is_default_thresholds(h.parse_thresholds([0.1, 0.2, 0.30000000000000004, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]))
	assert not h.is_default_thresholds(h.parse_thresholds(SMALL))
	assert h.segm_id_layout(h.DEFAULT_THRESHOLDS) == (1, 5)
	assert h.segm_id_layout(h.parse_thresholds(SMALL)) == (2, 7)
	# a fixed threshold: its position in the sorted list
	assert [h.fixed_level({'THRESHOLD_METHOD': t}) for t in (0.1, 0.5, 0.9, 'f', 'crossing')] == [1, 5, 9, None, None]
	assert h.fixed_level({'SEGM_THRESHOLDS': SMALL, 'THRESHOLD_METHOD': 0.06}) == 3
	assert h.fixed_level({'SEGM_THRESHOLDS': SMALL, 'THRESHOLD_METHOD': 0.1 + 1e-12}) == 5
	for t in (0.05, 0.5, 0.061):
		msg = rejects(h.fixed_level, {'SEGM_THRESHOLDS': SMALL, 'THRESHOLD_METHOD': t})
		assert '0.02, 0.04, 0.06, 0.08, 0.1' in msg, msg
	rejects(h.fixed_level, {'THRESHOLD_METHOD': 0.55})  # not one of the nine of 0.1-0.9
	assert h.thresholds_error({'SEGM_THRESHOLDS': SMALL, 'THRESHOLD_METHOD': 'f'}) is None
	assert h.thresholds_error({'SEGM_THRESHOLDS': [0, 0.5]}) is not None
	print('test_levels: OK')


def test_options():
	"""options_error: THRESHOLD_METHOD and SMALL_SEGM_MERGE checked before any step runs; segm_id_error."""
	ok = {'THRESHOLD_METHOD': 'f'}
	pal = h.PALIST
	h.PALIST = os.path.join(HERE, 'no_palist.csv')  # without palist.csv: the IDs are not checked
	try:
		assert h.options_error(ok) is None and h.options_error({}) is None
		for t in (True, '0.05', 'F', None, [0.5]):
			assert 'THRESHOLD_METHOD' in h.options_error({'THRESHOLD_METHOD': t}), t
		assert 'not one of the thresholds' in h.options_error({'THRESHOLD_METHOD': 0.05})
		assert h.options_error({'SEGM_THRESHOLDS': SMALL, 'THRESHOLD_METHOD': 0.06}) is None
		for m in ('similar', 'frequent', 'Similar', ' FREQUENT ', None, ''):
			assert h.options_error(dict(ok, SMALL_SEGM_MERGE=m)) is None, m
		for m in ('nearest', 1, ['similar']):
			assert 'SMALL_SEGM_MERGE' in h.options_error(dict(ok, SMALL_SEGM_MERGE=m)), m
		assert [h.merge_rule({'SMALL_SEGM_MERGE': m}) for m in ('Similar', None, '', 'frequent')] == \
			['similar', 'frequent', 'frequent', 'frequent']
		assert h.merge_rule({}) == 'frequent'
		with tempfile.TemporaryDirectory() as d:
			h.PALIST = os.path.join(d, 'palist.csv')
			open(h.PALIST, 'w').write('1\n1234567890\n')
			assert h.options_error(ok) is None  # 13-digit IDs with the default thresholds: 1234567890 9 99999
			msg = h.options_error({'SEGM_THRESHOLDS': SMALL})
			assert msg is None, msg  # 1234567890 05 9999999 = 1.2e18 < 2^63
			open(h.PALIST, 'w').write('9300000000\n')
			assert '9300000000059999999' in h.options_error({'SEGM_THRESHOLDS': SMALL})
	finally:
		h.PALIST = pal
	# the largest segm_id: area ID, the last k (or that of a fixed threshold) and the largest segment number
	assert h.segm_id_error('1120008380', SMALL) is None  # 10-digit HydroBASINS ID: 19 digits, fits int64
	assert h.segm_id_error('9223372036', SMALL, 1) is None  # 9223372036 01 9999999 < 2^63 - 1
	assert h.segm_id_error('9223372037', SMALL, 1) is not None  # 9223372037 01 9999999 > 2^63 - 1
	assert h.segm_id_error('9223372036853', h.DEFAULT_THRESHOLDS, 7) is None  # 9223372036853 7 99999 < 2^63 - 1
	assert h.segm_id_error('9223372036854', h.DEFAULT_THRESHOLDS, 7) is not None  # 9223372036854 7 99999 > 2^63 - 1
	assert h.segm_id_error('92233720368548', h.DEFAULT_THRESHOLDS) is not None  # 14 digits + 6
	assert h.segm_id_error('-1', SMALL) is None and h.segm_id_error('ES123', SMALL) is None
	assert abs(h.class_memory_gb(6000) - 2.52) < 0.01
	print('test_options: OK')


def test_set_option():
	"""The app writes NAME = value over the whole value, also over several lines, keeping the comment after it."""
	text = ('# comment\nSEGM_THRESHOLDS = [0.1, 0.2,\n                   0.3]  # inline\nSIMIL_CUTS = [10, 25]\n'
	        "SCENARIO = 'año'  # ñ\n")
	out = h.set_option(text, 'SEGM_THRESHOLDS', '[0.05, 0.1]')
	assert out == "# comment\nSEGM_THRESHOLDS = [0.05, 0.1]  # inline\nSIMIL_CUTS = [10, 25]\nSCENARIO = 'año'  # ñ\n", out
	out = h.set_option(text, 'SCENARIO', "'b'")
	assert out.endswith("SCENARIO = 'b'  # ñ\n") and 'SEGM_THRESHOLDS = [0.1, 0.2,\n' in out, out
	assert h.set_option(text, 'NEW', '1').endswith("# ñ\nNEW = 1\n")
	# a config.py that does not parse: the line of the option, as before
	assert h.set_option('A = 1\nB = (\n', 'A', '2') == 'A = 2\nB = (\n'
	print('test_set_option: OK')


def test_values():
	assert [h.threshold_label(v) for v in (0.1, 0.3, 3 / 10, 0.1 + 0.2, 0.02, 0.125, 1e-5)] == \
		['0.1', '0.3', '0.3', '0.3', '0.02', '0.125', '0.00001']
	with tempfile.TemporaryDirectory() as d:
		p = os.path.join(d, 'segm_thresholds.csv')
		h.write_run_thresholds(h.parse_thresholds(SMALL), p)
		assert open(p).read() == 'k value\n1 0.02\n2 0.04\n3 0.06\n4 0.08\n5 0.1\n'
		assert h.run_thresholds(p) == {1: 0.02, 2: 0.04, 3: 0.06, 4: 0.08, 5: 0.1}
		# the nine of 0.1-0.9 written give the same values as the fallback of the segmentations without the file
		h.write_run_thresholds(h.DEFAULT_THRESHOLDS, p)
		fallback = h.run_thresholds(os.path.join(d, 'missing.csv'), {})
		assert h.run_thresholds(p) == fallback == {k: k / 10 for k in range(1, 10)}
		# without the file, or with a row that cannot be read, the thresholds of config.py (not always k / 10)
		assert h.run_thresholds(os.path.join(d, 'missing.csv'), {'SEGM_THRESHOLDS': SMALL}) == dict(enumerate(SMALL, 1))
		open(p, 'w').write('k value\n1 0.02 x\n')
		assert h.run_thresholds(p, {'SEGM_THRESHOLDS': SMALL}) == dict(enumerate(SMALL, 1))
		assert h.run_thresholds(p, {'SEGM_THRESHOLDS': [0.5, 0.1]}) == fallback  # invalid option: the default ones
	print('test_values: OK')


def test_recorded():
	"""Default values of SEGM_THRESHOLDS and SMALL_SEGM_MERGE: not in the inputs of the segmentation."""
	base = {'THRESHOLD_METHOD': 'f'}
	segm = lambda **kw: {k: v for k, v in h.step_inputs('segm', dict(base, **kw), {}).items()
	                     if k in ('SEGM_THRESHOLDS', 'SMALL_SEGM_MERGE', 'thresholds segmented')}
	assert segm() == {}
	assert segm(SEGM_THRESHOLDS=list(h.DEFAULT_THRESHOLDS), SMALL_SEGM_MERGE='frequent') == {}
	assert segm(SEGM_THRESHOLDS=tuple(h.DEFAULT_THRESHOLDS)) == {}
	assert segm(SMALL_SEGM_MERGE='Frequent') == segm(SMALL_SEGM_MERGE=None) == {}  # read as 'frequent'
	assert segm(SMALL_SEGM_MERGE='similar') == segm(SMALL_SEGM_MERGE='Similar') == {'SMALL_SEGM_MERGE': "'similar'"}
	assert segm(SMALL_SEGM_MERGE='nearest') == {'SMALL_SEGM_MERGE': "'nearest'"}  # invalid: recorded as it is
	assert segm(SEGM_THRESHOLDS=SMALL) == {'SEGM_THRESHOLDS': '[0.02, 0.04, 0.06, 0.08, 0.1]'}
	assert segm(SEGM_THRESHOLDS=SMALL, THRESHOLD_METHOD=0.06)['thresholds segmented'] == '1-3'
	assert segm(THRESHOLD_METHOD=0.5) == {'thresholds segmented': '1-5'}
	assert segm(THRESHOLD_METHOD=0.55) == {}  # invalid: not recorded here, the scripts stop saying why
	print('test_recorded: OK')


if __name__ == '__main__':
	test_parse()
	test_levels()
	test_options()
	test_set_option()
	test_values()
	test_recorded()
	print('All the tests passed')
