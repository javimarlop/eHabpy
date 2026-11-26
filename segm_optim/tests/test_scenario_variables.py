#### eHabitat+ HFT model: test of the variables of a scenario for the similarity (similarity.py)
#### License: CC BY-SA 3.0
"""
Checks that similarity.py takes the variables of a scenario (ENV_VARS of its processing_metadata.txt) when their NoData
is float('nan'), float('inf') or float('-inf'), as setup.py writes them, keeping it NaN or infinite; that it stops with
an error, instead of taking those of config.py (which can be those of another scenario), when that file is there but
its ENV_VARS cannot be read; that it takes those of config.py only for results without that file; and that setup.py
writes every NoData so that it can be read back.

    cd segm_optim && python -s tests/test_scenario_variables.py
"""

import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim

import setup  # noqa: E402
import similarity  # noqa: E402

# config.py of the session: the 8 variables of another scenario
CFG = {'ENV_VARS': {f'other_{i}': {'file': f'other/{i}.tif', 'nodata': -9999} for i in range(1, 9)}}

# processing_metadata.txt as moranvar_plots.R writes it (config.py without its comments), with 7 variables
META = """eHabitat+ processing metadata: options used to produce the files in this folder
Date: 2026-11-26 10:00:00
Study area: SierraNevada_2000m_EPSG32630
Scenario: borreguiles_7vars

## config.py
GRASSDB = '~/grassdata'
RESOLUTION = 10
ENV_VARS = {
{entries}
}
SIMIL_CUTS = [10, 25, 50]

## moranvar_plots.R
upper = 1
PCA of ID 1: the segmentation used PC1-PC3 of 7 axes, 81.2 % of the variance (the first 3 axes, SEGM_PCA_VARIANCE = 0)
"""


def entries(nodata):
	"""The lines of ENV_VARS of 7 variables, as setup.py writes them, with the NoData NODATA (text) of each."""
	nodata = (nodata * 7)[:7]
	return '\n'.join(f"    {',' if i else ''}'set_B_{i + 1}_var': {{'file': 'set_B/{i + 1}_var{i + 1}.tif', 'nodata': {nd}}}"
	                 for i, nd in enumerate(nodata))


def scenario(d, text):
	root = os.path.join(d, 'scenario')
	shutil.rmtree(root, ignore_errors=True)
	os.makedirs(root)
	if text is not None:
		with open(os.path.join(root, 'processing_metadata.txt'), 'w', encoding='utf-8') as f:
			f.write(text)
	return root


def fails(root, *words):
	"""The error of scenario_variables for ROOT, which must say every one of WORDS."""
	try:
		got = similarity.scenario_variables(root, CFG)
	except ValueError as e:
		msg = str(e)
		assert all(w in msg for w in words), (msg, words)
		return msg
	raise AssertionError(f'no error, took {got}')


def main():
	d = tempfile.mkdtemp()
	try:
		# NoData NaN (the case seen on 2026-11-26): the 7 variables of the scenario, not the 8 of config.py
		env, source = similarity.scenario_variables(scenario(d, META.replace('{entries}', entries(["float('nan')"]))), CFG)
		assert source == 'processing_metadata.txt', source
		assert list(env) == [f'set_B_{i}_var' for i in range(1, 8)], list(env)
		assert env['set_B_3_var']['file'] == 'set_B/3_var3.tif'
		assert all(isinstance(v['nodata'], float) and math.isnan(v['nodata']) for v in env.values()), env
		pairs = similarity.variable_pairs(env, os.path.join(d, 'nowhere'))
		assert [p['key'] for p in pairs] == [f'var{i}' for i in range(1, 8)], pairs

		# infinite NoData, numbers, None, spaces around = and a closing } with spaces before the next option
		text = META.replace('{entries}', entries(["float('inf')", "float('-inf')", "float( ' -Infinity ' )", '-9999',
		                                          '0.5', 'None', "float('NaN')"])).replace('ENV_VARS = {', 'ENV_VARS={')
		env, source = similarity.scenario_variables(scenario(d, text), CFG)
		nd = [v['nodata'] for v in env.values()]
		assert source == 'processing_metadata.txt' and nd[:6] == [math.inf, -math.inf, -math.inf, -9999, 0.5, None], nd
		assert math.isnan(nd[6]), nd

		# ENV_VARS that cannot be read: an error with the file and the line, never the variables of config.py
		root = scenario(d, META.replace('{entries}', entries(['-9999', "float('nan')", 'nan'])))  # nan: not Python
		msg = fails(root, 'processing_metadata.txt', 'line 12', 'config.py')  # ENV_VARS = { is line 9
		assert 'other_1' not in msg
		fails(scenario(d, META.replace('{entries}', entries(["__import__('os').getcwd()"]))), 'line 10')
		fails(scenario(d, META.replace('{entries}', entries(["float('1')"]))), 'not a number')  # only nan, inf, -inf
		fails(scenario(d, META.replace('{entries}', entries(['-9999']).replace("'nodata':", "'nodata'::", 1))),
		      'line 10')  # a syntax error
		fails(scenario(d, META.replace('ENV_VARS = {\n{entries}\n}', "ENV_VARS = {\n    'a', 'b'\n}")),
		      'not a dictionary')
		fails(scenario(d, META.replace('ENV_VARS = {\n{entries}\n}', '')), 'there are none')
		fails(scenario(d, META.replace('{entries}\n}', entries(['-9999']) + '}')), 'as setup.py writes them')

		# results without processing_metadata.txt (made before it): those of config.py
		env, source = similarity.scenario_variables(scenario(d, None), CFG)
		assert source == 'config.py' and env == CFG['ENV_VARS'], (source, env)

		# literal_value: the literals of ast.literal_eval and the non-finite floats only
		assert similarity.literal_value("[1, -2.5, 'a', None, (True,)]") == [1, -2.5, 'a', None, (True,)]
		assert similarity.literal_value("-float('inf')") == -math.inf
		for bad in ("float(x)", "float('nan', 1)", "float(s='nan')", "int('nan')", "x.float('nan')", "open('f')"):
			try:
				similarity.literal_value(bad)
			except ValueError:
				continue
			raise AssertionError(f'{bad} was evaluated')

		# setup.py writes every NoData so that config.py and processing_metadata.txt can be read back
		assert setup.fmt_num(math.inf) == "float('inf')" and setup.fmt_num(-math.inf) == "float('-inf')"
		for v in (math.nan, math.inf, -math.inf, -9999.0, 0.5, -3.4028234663852886e+38, None):
			back = similarity.literal_value(setup.fmt_num(v))
			assert back == v or (v != v and back != back), (v, setup.fmt_num(v), back)
	finally:
		shutil.rmtree(d)
	print('All the tests passed')


if __name__ == '__main__':
	main()
