#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 (2024)
#### License: CC BY-SA 3.0
####
#### Per-HFT SUMMARY (MEAN or MEDIAN according to config.py) + variance summary
#### of the input variables for each segmentation threshold level k (1..n, the thresholds of SEGM_THRESHOLDS;
#### 1..9 for the nine of 0.1-0.9).
#### ehab_segm_core.py computes them for all the segments of an area and threshold at once, from the
#### segment raster tiffs/segm_<id>_<k>.tif written by segmentation_pca_par.py (one NumPy pass per variable).
####
#### Usage:
####     python getmean_median_var.py
####     OR
####     from getmean_median_var import *
####     run_batch_all()

import sys
import numpy as np
import ehab_segm_core as _core
import hft_state
from ehab_segm_core import safelyMakeDir, initglobalmaps

# Importar la opción de agregación desde config.py
try:
    import config
    STAT_AGG = getattr(config, 'STAT_AGG', 'mean')
except ImportError:
    print("WARNING: No se encontró config.py. Se usará 'mean' por defecto.")
    STAT_AGG = 'mean'

# Seleccionar la función de agregación de NumPy segun config
if str(STAT_AGG).lower() == 'median':
    _AGG = np.median
    print("--> Modo de agregación seleccionado: MEDIANA (numpy.median)")
else:
    _AGG = np.mean
    print("--> Modo de agregación seleccionado: MEDIA (numpy.mean)")


def run_batch_all():
	# Stops saying what to run if the segmentation is out of date with config.py, palist.csv or the data
	hft_state.check_ready('stats')
	hft_state.start('stats')
	out = _core.run_batch_all(_AGG)
	hft_state.done('stats')
	return out


def run_park(paid):
	hft_state.check_ready('stats')
	return _core.run_park(paid, _AGG)


def run_batch(k):
	return _core.run_batch(k, _AGG)


def ehabitat(ecor, nw, nwpathout, k):
	return _core.ehabitat(ecor, nw, nwpathout, k, _AGG)


# --- Backward-compatible per-level wrappers (ehabitat1..9 / run_batch1..9, the nine thresholds of 0.1-0.9; any level
# k: ehabitat(ecor, nw, nwpathout, k) and run_batch(k)) ------
def _make_level(k):
	def _eh(ecor, nw, nwpathout, _k=k):
		return _core.ehabitat(ecor, nw, nwpathout, _k, _AGG)

	def _rb(_k=k):
		return _core.run_batch(_k, _AGG)
	return _eh, _rb


for _k in range(1, 10):
	globals()['ehabitat' + str(_k)], globals()['run_batch' + str(_k)] = _make_level(_k)
del _k


if __name__ == '__main__':
	run_batch_all()
