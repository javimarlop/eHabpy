"""Folders of the copy of eHabpy the app works on (the one it belongs to)."""

import json
import os
import sys
import threading

APP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # app/
ROOT = os.path.dirname(APP_DIR)  # the copy of eHabpy
SEGM = os.path.join(ROOT, 'segm_optim')
INVARS = os.path.join(ROOT, 'inVars')
PA_DIR = os.path.join(SEGM, 'pa')
RESULTS = os.path.join(SEGM, 'results')
UI = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ui')
SETTINGS = os.path.join(SEGM, 'app_settings.json')  # choices of the app that config.py does not keep
APP_RECIPE = os.path.join(SEGM, 'app_recipe.yaml')  # last recipe of MonitorEO run from the app

ECO_FOLDER = 'ecoregion'  # subfolder of a set of variables with those of the ecoregion scale (as get_vars_gee.py)
EXCLUDED = '_excluded'  # subfolder of the variables left out of a set (setup.py only reads the files of the folder)

# the scripts of eHabpy (hft_state, setup, get_vars_gee) are imported from segm_optim
if SEGM not in sys.path:
	sys.path.insert(0, SEGM)

_lock = threading.Lock()


def load_settings():
	try:
		with open(SETTINGS, encoding='utf-8') as f:
			return json.load(f)
	except (OSError, ValueError):
		return {}


def save_settings(**changes):
	"""Update the settings of the app with CHANGES (None removes a key) and return them."""
	with _lock:
		s = load_settings()
		for k, v in changes.items():
			if v is None:
				s.pop(k, None)
			else:
				s[k] = v
		tmp = SETTINGS + '.tmp'
		with open(tmp, 'w', encoding='utf-8') as f:
			json.dump(s, f, indent=1, ensure_ascii=False)
		os.replace(tmp, SETTINGS)
		return s


def rel(path):
	"""PATH relative to the copy of eHabpy when it is inside it ('/' as separator), else PATH."""
	try:
		r = os.path.relpath(path, ROOT)
	except ValueError:  # another drive on Windows
		return path
	return path if r.startswith('..') else r.replace(os.sep, '/')


def inside(path, folder):
	"""Whether PATH is FOLDER or inside it."""
	try:
		return os.path.commonpath([os.path.realpath(path), os.path.realpath(folder)]) == os.path.realpath(folder)
	except ValueError:
		return False
