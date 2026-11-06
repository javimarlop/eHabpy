#### Author: Javier Martinez-Lopez (UTF-8) 2014 - 2021
#### Modernized for Python 3 (2024)
#### License: CC BY-SA 3.0
####
#### Per-HFT MEDIAN + variance summary of the 9 input variables, for each of the
#### nine segmentation threshold levels (1..9).
####
#### Same as getmeanvar.py but the per-variable summary uses ``numpy.median``
#### instead of ``numpy.mean`` (the variance columns still use ``numpy.var``).
#### The shared logic lives in ``ehab_segm_core``.
####
#### Usage (unchanged):
####     python
####     from getmedianvar import *
####     run_batch_all()
####     exit()

import numpy as np
import ehab_segm_core as _core
from ehab_segm_core import safelyMakeDir, initglobalmaps

_AGG = np.median


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
