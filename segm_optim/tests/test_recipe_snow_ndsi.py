#### eHabitat+ desktop app: test of the snow by its NDSI in the recipes of MonitorEO (quality snow_ndsi)
#### License: CC BY-SA 3.0
"""
Checks that the recipe of the form keeps quality snow_ndsi (MonitorEO 0.3.80 or later): written when it is a number,
left out when it is off (so that a version of monitoreo without it does not refuse the recipe), and that the recipe
of the repository does not have it on. Without Earth Engine nor network.

    cd segm_optim && python -s tests/test_recipe_snow_ndsi.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), 'app'))  # ehabitat_app

from ehabitat_app import gee  # noqa: E402


def test_snow_ndsi():
	rec = gee.default_recipe()
	assert rec['quality']['snow_ndsi'] is None  # the form shows it off
	assert 'snow_ndsi' not in gee.clean_recipe(rec)['quality']
	rec['quality']['snow_ndsi'] = ''
	assert 'snow_ndsi' not in gee.clean_recipe(rec)['quality']
	rec['quality']['snow_ndsi'] = 0.4
	assert gee.clean_recipe(rec)['quality']['snow_ndsi'] == 0.4
	rec['quality']['topographic'] = True  # an option the form does not show stays
	assert gee.clean_recipe(rec)['quality']['topographic'] is True
	print('test_snow_ndsi: OK')


if __name__ == '__main__':
	test_snow_ndsi()
	print('All the tests passed')
