#### eHabitat+: test of the warning on the buffer of the study grid of a download (get_vars_gee.buffer_warning)
#### License: CC BY-SA 3.0
"""
Checks that a download warns when the buffer of its study grid (study: buffer in the recipe) is smaller than that of
the model (PA_BUFFER of config.py, "Buffer around the study area" in HFT mapping), and not when it is not given (the
download takes PA_BUFFER), equal or larger; and that the app reads PA_BUFFER as the model does (0 when it is not in
config.py). Without Earth Engine.

    cd segm_optim && python -s tests/test_buffer_warning.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'app'))  # ehabitat_app

import get_vars_gee as g  # noqa: E402


def test_buffer_warning():
	w = g.buffer_warning({'study': {'resolution': 10, 'buffer': 100}}, 500)
	assert w.startswith('WARNING: the buffer of the study grid (100 m) is smaller than that of the model (PA_BUFFER = '
	                    '500 m'), w
	for rec, pa in (({'study': {'buffer': None}}, 500), ({'study': {}}, 500), ({}, 500), ({'study': {'buffer': 500}}, 500),
	                ({'study': {'buffer': 600}}, 500), ({'study': {'buffer': 0}}, 0), ({'study': {'buffer': 'x'}}, 500),
	                ({'study': {'buffer': 100}}, None), ({'study': {'buffer': True}}, 500)):
		assert g.buffer_warning(rec, pa) is None, (rec, pa)
	print('test_buffer_warning: OK')


def test_app_pa_buffer():
	from ehabitat_app import model
	real = model.read_config
	try:
		for cfg, want in (({'PA_BUFFER': 500}, 500.0), ({'PA_BUFFER': 0}, 0.0), ({}, 0.0), ({'PA_BUFFER': 'x'}, 0.0),
		                  ({'PA_BUFFER': 12.5}, 12.5)):
			model.read_config = lambda cfg=cfg: cfg
			assert model.pa_buffer() == want, (cfg, model.pa_buffer())
	finally:
		model.read_config = real
	print('test_app_pa_buffer: OK')


if __name__ == '__main__':
	test_buffer_warning()
	test_app_pa_buffer()
	print('All the tests passed')
