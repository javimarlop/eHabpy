#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 (2024)
#### License: CC BY-SA 3.0
####
#### Per-HFT SUMMARY (MEAN or MEDIAN according to config.py) + variance summary 
#### of the input variables for each of the nine segmentation threshold levels (1..9).
####
#### Usage:
####     python getmean_median_var.py
####     OR
####     from getmean_median_var import *
####     run_batch_all()

import sys
import numpy as np
import ehab_segm_core as _core
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
	return _core.run_batch_all(_AGG)


def run_park(paid):
	return _core.run_park(paid, _AGG)


def run_batch(k):
	return _core.run_batch(k, _AGG)


def ehabitat(ecor, nw, nwpathout, k):
	return _core.ehabitat(ecor, nw, nwpathout, k, _AGG)


# --- Backward-compatible per-level wrappers (ehabitat1..9 / run_batch1..9) ------
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
