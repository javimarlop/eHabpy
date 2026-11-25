#### eHabitat+: the order of the variables of a set
#### License: CC BY-SA 3.0

"""
The order of the variables of a set (a folder of ../inVars) is that of the model: the first one defines the pixels
analysed (FIRST_VAR of segmentation_pca_par.py) and they are the inputs of the PCA in that order. setup.py writes it in
config.py (ENV_VARS0, ENV_VARS) as:

  1. the files listed in ORDER_FILE of the folder (one name per line; written by the app when the variables are moved
     up or down), in that order;
  2. then those already in config.py (ENV_VARS: the order of the last setup.py), in that order;
  3. then the others, in the order of their names (sorted(), so capitals before small letters).

The files keep their names: renaming them would make the GRASS maps (named after them) be imported again, and a new
download of the set (gee_vars.json) download them again with their first names.
"""

import os

ORDER_FILE = 'variables_order.txt'
ORDER_CHANGED = 'the order of the variables changed'  # the reason of setup_changes (hft_state.py)


def only_order(reasons):
	"""Whether the REASONS of setup_changes to run setup.py are only a new order of the variables (it imports nothing
	then: the maps are in GRASS already, it only writes config.py)."""
	return bool(reasons) and all(r.startswith(ORDER_CHANGED) for r in reasons)


def read_order(folder):
	"""The names of the files listed in ORDER_FILE of FOLDER, in their order ([] without it)."""
	try:
		with open(os.path.join(folder, ORDER_FILE), encoding='utf-8-sig') as f:  # a BOM of Notepad is left out
			lines = [x.strip() for x in f]
	except OSError:
		return []
	out = []
	for x in lines:
		if x and not x.startswith('#') and os.path.basename(x) == x and x not in out:
			out.append(x)
	return out


def write_order(folder, names):
	"""Write the order of the files NAMES (basenames) to ORDER_FILE of FOLDER."""
	text = ('# The order of the variables of this set for eHabitat+ (the first one defines the pixels analysed): one\n'
	        '# file a line. The files not listed go after these, as setup.py found them (see segm_optim/var_order.py).\n'
	        + ''.join(f'{n}\n' for n in names))
	path = os.path.join(folder, ORDER_FILE)
	with open(path + '.tmp', 'w', encoding='utf-8') as f:
		f.write(text)
	os.replace(path + '.tmp', path)


def ordered(paths, folder, previous=None):
	"""PATHS (files of FOLDER) in the order of the model: those of ORDER_FILE, then those of PREVIOUS (basenames in the
	order of the last setup.py), then the others sorted by name."""
	by_name = {os.path.basename(p): p for p in paths}
	out = []
	for n in read_order(folder) + list(previous or []):
		if n in by_name and by_name[n] not in out:
			out.append(by_name[n])
	return out + sorted(p for p in paths if p not in out)


def previous_order(env_vars, folder):
	"""The basenames of the files of the set FOLDER (relative to ../inVars, '' for inVars itself) in ENV_VARS of
	config.py, in their order."""
	out = []
	for v in (env_vars or {}).values():
		f = str((v or {}).get('file') or '').replace('\\', '/')
		d, name = (f.rsplit('/', 1) if '/' in f else ('', f))
		if d == folder and name:
			out.append(name)
	return out
