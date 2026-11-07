"""Folders of the copy of eHabpy the app works on (the one it belongs to)."""

import json
import os
import shutil
import subprocess
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


def wsl():
	"""Whether this is Linux in WSL (the Windows Subsystem for Linux), where the web browser and the File Explorer are
	those of Windows."""
	if not sys.platform.startswith('linux'):
		return False
	if os.environ.get('WSL_DISTRO_NAME') or os.environ.get('WSL_INTEROP'):
		return True
	try:
		with open('/proc/sys/kernel/osrelease', encoding='utf-8') as f:
			return 'microsoft' in f.read().lower()
	except OSError:
		return False


def open_in_windows(target):
	"""In WSL, open TARGET with Windows: a web address in its web browser or a folder of Linux in its File Explorer
	(as \\\\wsl.localhost\\...). Returns False when Windows cannot be reached from Linux."""
	exe = shutil.which('explorer.exe') or '/mnt/c/Windows/explorer.exe'
	if not os.path.exists(exe):
		return False
	try:
		if not target.startswith(('http://', 'https://')):
			r = subprocess.run(['wslpath', '-w', target], capture_output=True, text=True)
			if r.returncode != 0 or not r.stdout.strip():
				return False
			target = r.stdout.strip()
		# explorer.exe ends with exit code 1 even when it opened TARGET: not checked
		subprocess.Popen([exe, target], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	except OSError:
		return False
	return True


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
