#### eHabitat+: test of the plan of a download from MonitorEO (get_vars_gee.py --dry-run)
#### License: CC BY-SA 3.0
"""
Checks that the plan of get_vars_gee.py --dry-run gives the number of images of each variable with its unit and period
("2,027 images (2018-01-01 to 2021-12-31)", "fixed image"), warns when a variable has no images or very few, and ends
with the size of the download by grid and in all (files × MB per variable). Without Earth Engine: the connection and
the count of the images are replaced, and the grids are those of monitoreo.grid.Grid.

    cd segm_optim && python -s tests/test_dry_run_plan.py
"""

import contextlib
import io
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # segm_optim

import get_vars_gee as g  # noqa: E402
from monitoreo import recipe  # noqa: E402
from monitoreo.grid import Grid  # noqa: E402
from osgeo import ogr, osr  # noqa: E402


class Out(io.StringIO):
	def reconfigure(self, **kw):  # as sys.stdout (main sets line_buffering)
		pass


class FakePlan:
	def __init__(self, specs, grids):
		self.specs, self.grids = specs, grids


def test_texts():
	spec = {'file': '1_NDVI_S2_p95.tif', 'period': ['2018-01-01', '2021-12-31'], 'season': None}
	assert g.images_text(spec, 2027) == '2,027 images (2018-01-01 to 2021-12-31)'
	assert g.images_text(spec, 1, True).startswith('1 image (') and 'WARNING: very few images' in g.images_text(spec, 1, True)
	assert 'WARNING' not in g.images_text(spec, 4)  # a collection of the user: maybe one image a year
	assert g.images_text(spec, None) == 'fixed image'
	assert g.images_text(spec, 'ERROR: no data') == 'ERROR: no data'
	summer = dict(spec, season=['06-01', '09-30'])
	t = g.images_text(summer, 0)
	assert t.startswith('0 images (2018-01-01 to 2021-12-31, season 06-01 to 09-30)') and \
		'WARNING: no images in the period and season' in t, t
	assert 'WARNING' not in g.images_text(spec, g.FEW_IMAGES)
	assert g.mb_text(0.5) == '0.50 MB' and g.mb_text(64.4) == '64 MB' and g.mb_text(2048) == '2.0 GB'
	assert g.mb_text(1023.6) == '1.0 GB' and g.mb_text(0.999) == '1 MB' and g.mb_text(0.001) == '0.01 MB'
	assert g.mb_text(3.3 * 1024 ** 2) == '3.3 TB' and g.mb_text(0) == '0 MB'

	# 6 variables on two grids (64 MB and 83 MB per variable): about 0.9 GB in 12 files
	study = Grid('EPSG:32630', 10, 440000, 4120000, 4096, 4096)  # 16 777 216 pixels: 64 MB
	eco = Grid('EPSG:32630', 100, 300000, 4300000, 4096, 5312)  # 21 757 952 pixels: 83 MB
	plan = FakePlan([{'file': f'{i}_v.tif'} for i in range(1, 7)], {'study': study, 'ecoregion': eco})
	lines = g.total_text(plan)
	assert lines[0].startswith('Size of all the files of the plan'), lines
	assert lines[1] == '  study grid: 6 files, about 384 MB', lines
	assert lines[2] == '  ecoregion grid: 6 files, about 498 MB', lines
	assert lines[-1] == '  in all: about 882 MB in 12 files', lines
	# with the study area at the pixel of the ecoregion grid and the true colour (three bands on the study grid)
	small = Grid('EPSG:32630', 100, 440000, 4120000, 410, 410)
	lines = g.total_text(plan, small, {'index': 'NDVI'})
	assert '  study area at the pixel of the ecoregion grid: 6 files, about 4 MB' in lines, lines
	assert '  true colour: 1 file, about 192 MB' in lines, lines
	assert lines[-1] == '  in all: about 1.1 GB in 19 files', lines
	assert g.total_text(plan, area=True)[2] == '  comparison area grid: 6 files, about 498 MB'
	print('test_texts: OK')


def test_dry_run():
	"""get_vars_gee.py --dry-run with the recipe of the repository: its last lines."""
	tmp = tempfile.mkdtemp()
	real = (g.update_monitoreo, recipe.n_images, recipe.plan, sys.argv)
	try:
		srs = osr.SpatialReference()
		srs.ImportFromEPSG(32630)
		ds = ogr.GetDriverByName('GPKG').CreateDataSource(os.path.join(tmp, 'area.gpkg'))
		lyr = ds.CreateLayer('area', srs, ogr.wkbPolygon)
		f = ogr.Feature(lyr.GetLayerDefn())
		f.SetGeometry(ogr.CreateGeometryFromWkt('POLYGON((450000 4100000,460000 4100000,460000 4110000,450000 4100000))'))
		lyr.CreateFeature(f)
		ds = None
		import monitoreo.services as services
		real_init = services.init
		services.init = lambda project=None: None

		def plan(rec, area, layer=None, buffer=None, only=None, extra=(), warn=print):
			rec = recipe.read(rec, extra)
			study = Grid.around('EPSG:32630', float(rec['study']['resolution']), 450000, 4100000, 460000, 4110000,
			                    margin=(buffer or 0) + 20)
			return recipe.Plan(rec, recipe.expand(rec, warn), {'study': study}, os.path.basename(area))

		recipe.plan = plan
		g.update_monitoreo = lambda check=True: None
		counts = iter([2027, 0, 3])
		recipe.n_images = lambda p: {s['file']: (None if recipe.variables.kind(s) != 'catalog' else next(counts)) for s in p.specs}
		sys.argv = ['get_vars_gee.py', '--dry-run', '--no-update', '--project', 'none', '--only', 'study',
		            '--study-area', os.path.join(tmp, 'area.gpkg'), '--recipe', os.path.join(os.path.dirname(HERE),
		                                                                                 'gee_vars.yaml')]
		out = Out()
		with contextlib.redirect_stdout(out):
			g.main()
		text = out.getvalue()
		lines = text.splitlines()
		i = lines.index('Images of each variable over the study area:')
		images = lines[i + 1:]
		assert images[0].endswith(': 2,027 images (2018-01-01 to 2021-12-31)'), images
		assert 'WARNING: no images in the period' in images[1], images
		assert 'WARNING: very few images' in images[2], images  # a variable of the catalog
		assert all(x.endswith(': fixed image') for x in images[3:6]), images  # slope, temp, precip
		n = sum(1 for x in images if x.startswith('  ') and '.tif: ' in x)
		j = lines.index('Size of all the files of the plan (uncompressed; those already up to date are not downloaded '
		                'again):')
		assert lines[j + 1].startswith(f'  study grid: {n} files, about '), lines[j:]
		assert lines[-1].startswith('  in all: about ') and lines[-1].endswith(f' in {n} files'), lines[-1]
		print('test_dry_run: OK')
		print('\n'.join(lines[i:]))
	finally:
		g.update_monitoreo, recipe.n_images, recipe.plan, sys.argv = real
		try:
			services.init = real_init
		except NameError:
			pass
		shutil.rmtree(tmp)


if __name__ == '__main__':
	test_texts()
	test_dry_run()
	print('All the tests passed')
