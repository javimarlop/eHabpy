"""Activation of the conda environment of eHabitat+ when the app is not opened from an activated one (the shortcut of
the Desktop runs its Python directly): the scripts of the model need GRASS GIS, Rscript and the GDAL tools in the PATH,
and the variables that the activation scripts of conda set (GDAL_DATA, PROJ_DATA, ...).

The activated environment is asked once to micromamba (or conda) with 'run' and kept in the environment
(.ehabitat_activation.json) until a package changes; without micromamba or conda, the folders of the environment are
added to the PATH."""

import json
import os
import shutil
import subprocess
import sys

CACHE = '.ehabitat_activation.json'
DUMP = 'import json, os; print(json.dumps(dict(os.environ)))'
SKIP = {'PWD', 'OLDPWD', '_', 'SHLVL', 'PS1', 'PROMPT', 'TERM_SESSION_ID', 'CONDA_PROMPT_MODIFIER'}
NO_WINDOW = getattr(subprocess, 'CREATE_NO_WINDOW', 0)


def _same(a, b):
	try:
		return os.path.normcase(os.path.realpath(a)) == os.path.normcase(os.path.realpath(b))
	except (OSError, ValueError):
		return False


def bin_dirs(prefix=sys.prefix):
	"""Folders of the programs of a conda environment, as conda puts them in the PATH."""
	if os.name == 'nt':
		return [prefix] + [os.path.join(prefix, *p) for p in (('Library', 'mingw-w64', 'bin'), ('Library', 'usr', 'bin'),
		                                                       ('Library', 'bin'), ('Scripts',), ('bin',))]
	return [os.path.join(prefix, 'bin')]


def _mamba_candidates(prefix):
	"""Executables of micromamba or conda that can run the environment PREFIX."""
	exe = '.exe' if os.name == 'nt' else ''
	found = [os.environ.get('EHAB_MAMBA'), os.environ.get('MAMBA_EXE'),
	         os.path.join(os.path.dirname(prefix), 'micromamba', 'micromamba' + exe),  # the installer of the app
	         os.environ.get('CONDA_EXE')]
	found += [shutil.which(n) for n in ('micromamba', 'mamba', 'conda')]
	out = []
	for f in found:
		if f and os.path.isfile(f) and f not in out:
			out.append(f)
	return out


def _history_sig(prefix):
	try:
		st = os.stat(os.path.join(prefix, 'conda-meta', 'history'))
		return f'{st.st_size} {st.st_mtime_ns}'
	except OSError:
		return None


def _ask(exe, prefix):
	"""Environment variables of PREFIX activated by EXE (micromamba or conda), or None."""
	env = dict(os.environ)
	if os.path.basename(exe).lower().startswith('micromamba'):
		env.setdefault('MAMBA_ROOT_PREFIX', os.path.dirname(exe))
	py = 'python.exe' if os.name == 'nt' else 'python'
	try:
		r = subprocess.run([exe, 'run', '-p', prefix, py, '-c', DUMP], capture_output=True, text=True, timeout=120,
		                   env=env, stdin=subprocess.DEVNULL, creationflags=NO_WINDOW)
	except (OSError, subprocess.SubprocessError):
		return None
	if r.returncode != 0:
		return None
	for line in reversed(r.stdout.strip().splitlines()):
		try:
			found = json.loads(line)
		except ValueError:
			continue
		if isinstance(found, dict) and found.get('PATH'):
			return found
	return None


def activate(prefix=sys.prefix):
	"""Make the processes started by the app see the environment PREFIX as activated. Returns how: 'active' (it was
	already), 'cache', 'micromamba' or 'conda' (asked to them) or 'path' (only its folders in the PATH)."""
	os.environ.setdefault('PYTHONNOUSERSITE', '1')
	os.environ['PYTHONUNBUFFERED'] = '1'  # the scripts write their messages to the log of the app as they go
	os.environ.setdefault('PYTHONIOENCODING', 'utf-8')
	if os.environ.get('CONDA_PREFIX') and _same(os.environ['CONDA_PREFIX'], prefix):
		return 'active'
	how, found = None, None
	cache = os.path.join(prefix, CACHE)
	sig = _history_sig(prefix)
	try:
		with open(cache, encoding='utf-8') as f:
			c = json.load(f)
		if sig and c.get('history') == sig:
			found, how = c['added'], 'cache'
	except (OSError, ValueError, KeyError):
		pass
	if found is None:
		for exe in _mamba_candidates(prefix):
			env = _ask(exe, prefix)
			if env:
				# only what the activation adds or changes (the rest is that of this computer, at this moment)
				found = {k: v for k, v in env.items() if k not in SKIP and os.environ.get(k) != v}
				how = 'micromamba' if 'micromamba' in os.path.basename(exe).lower() else 'conda'
				try:
					with open(cache, 'w', encoding='utf-8') as f:
						json.dump({'history': sig, 'added': found}, f)
				except OSError:
					pass
				break
	if found is not None:
		path = found.pop('PATH', None)
		os.environ.update(found)
		if path:
			# the folders of the environment first, then those of this computer now
			keep = [p for p in os.environ.get('PATH', '').split(os.pathsep) if p and p not in path.split(os.pathsep)]
			os.environ['PATH'] = os.pathsep.join([path] + keep)
	else:
		how = 'path'
	# in every case, the folders of the environment at the start of the PATH
	dirs = [d for d in bin_dirs(prefix) if os.path.isdir(d)]
	rest = [p for p in os.environ.get('PATH', '').split(os.pathsep) if p and p not in dirs]
	os.environ['PATH'] = os.pathsep.join(dirs + rest)
	os.environ.setdefault('CONDA_PREFIX', prefix)
	if how == 'path':
		share = os.path.join(prefix, 'Library', 'share') if os.name == 'nt' else os.path.join(prefix, 'share')
		for var, sub in (('GDAL_DATA', 'gdal'), ('PROJ_DATA', 'proj'), ('PROJ_LIB', 'proj')):
			if os.path.isdir(os.path.join(share, sub)):
				os.environ.setdefault(var, os.path.join(share, sub))
	return how


def tools():
	"""Versions of the programs the model needs: {name: version or None}."""
	def first_line(cmd, shell=False):
		try:
			r = subprocess.run(cmd, capture_output=True, text=True, timeout=60, shell=shell,
			                   stdin=subprocess.DEVNULL, creationflags=NO_WINDOW)
		except (OSError, subprocess.SubprocessError):
			return None
		out = (r.stdout or r.stderr or '').strip().splitlines()
		return out[0].strip() if r.returncode == 0 and out else None

	out = {'python': sys.version.split()[0]}
	try:
		from osgeo import gdal
		out['gdal'] = gdal.__version__
	except ImportError:
		out['gdal'] = None
	grass = None
	for exe in ('grass', 'grass85', 'grass84', 'grass83'):
		v = first_line([exe, '--version'], shell=(os.name == 'nt'))
		if v:
			grass = v.replace('GRASS GIS', '').strip() or v
			break
	out['grass'] = grass
	r = first_line(['Rscript', '-e', 'cat(R.version$major, R.version$minor, sep=".")'])
	out['r'] = r
	try:
		import importlib.metadata
		out['monitoreo'] = importlib.metadata.version('monitoreo')
	except Exception:
		out['monitoreo'] = None
	try:
		import webview
		out['pywebview'] = getattr(webview, '__version__', None) or importlib.metadata.version('pywebview')
	except Exception:
		out['pywebview'] = None
	return out
