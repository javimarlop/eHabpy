"""Updates of eHabitat+ from its GitHub repository, so that an installed copy always has the last version.

Every time the app opens (unless EHAB_NO_UPDATE=1 or auto_update is false in its settings), before its window, it
compares its copy with the last commit of its branch in GitHub and, if there is a newer one, updates it:

- a copy without git (downloaded as a zip by the installer or from GitHub) gets the files of the new commit from the
  zip of that commit. The files that the user changed since they were installed (config.py, gee_vars.yaml, ...) are
  kept, and so is everything that is not a file of the repository (inVars, results, ...): the hash of every file
  installed is kept in .ehabitat_version.json, with the commit. The files deleted in the repository are deleted if
  the user did not change them. A zip without the app (app/ehabitat.py) is never applied, nor an older commit;
- a git clone is updated with 'git pull --ff-only' of the branch it tracks: git never overwrites local changes, and a
  clone with commits of its own or changes in the same files is left as it is (the app says so).

After an update the app starts again with the new code. When environment.yml changed, the conda environment is
updated in the background (micromamba install, pip for its pip packages and, on Windows, CRAN for the R packages
marked [win: CRAN name] in it, which conda-forge does not have for Windows) and the app asks to restart when done.
"""

import hashlib
import importlib.metadata
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
import urllib.request
import zipfile

from . import envsetup, jobs, paths

REPO = 'javimarlop/eHabpy'
API = f'https://api.github.com/repos/{REPO}'
RECORD = os.path.join(paths.ROOT, '.ehabitat_version.json')
REQUIRED = 'app/ehabitat.py'  # a version without the app is never applied
# files of the repository that the user edits: kept when there is no record of what was installed
PROTECTED = {'segm_optim/config.py', 'segm_optim/gee_vars.yaml', 'segm_optim/palist.csv'}
NEVER = ('.git/', '.ehabitat_version.json')
NO_WINDOW = getattr(subprocess, 'CREATE_NO_WINDOW', 0)

_state = {'checked': None, 'result': None, 'restart': False, 'env_job': None, 'error': None}
_lock = threading.Lock()


def _say(log, text):
	(log or (lambda t: print(t, flush=True)))(text)


def _get(url, accept='application/vnd.github+json', timeout=10):
	req = urllib.request.Request(url, headers={'Accept': accept, 'User-Agent': 'eHabitat'})
	with urllib.request.urlopen(req, timeout=timeout) as r:
		return r.read()


def _sha1(data):
	return hashlib.sha1(data).hexdigest()


def _file_sha1(path):
	with open(path, 'rb') as f:
		return _sha1(f.read())


def load_record():
	try:
		with open(RECORD, encoding='utf-8') as f:
			return json.load(f)
	except (OSError, ValueError):
		return {}


def _save_record(rec):
	tmp = RECORD + '.tmp'
	with open(tmp, 'w', encoding='utf-8') as f:
		json.dump(rec, f, indent=1, sort_keys=True)
	os.replace(tmp, RECORD)


def mode():
	return 'git' if os.path.exists(os.path.join(paths.ROOT, '.git')) else 'archive'


def _git(*args, timeout=60):
	env = dict(os.environ, GIT_TERMINAL_PROMPT='0', GIT_SSH_COMMAND='ssh -o BatchMode=yes')
	try:
		r = subprocess.run(['git', '-C', paths.ROOT, *args], capture_output=True, text=True, timeout=timeout,
		                   env=env, stdin=subprocess.DEVNULL, creationflags=NO_WINDOW)
	except (OSError, subprocess.SubprocessError) as e:
		return None, str(e)
	return (r.stdout.strip() if r.returncode == 0 else None), r.stderr.strip()


def branch():
	"""Branch of GitHub followed by a copy without git: that of its record, EHAB_BRANCH or master."""
	return load_record().get('branch') or os.environ.get('EHAB_BRANCH') or 'master'


def enabled():
	if os.environ.get('EHAB_NO_UPDATE', '') not in ('', '0'):
		return False
	return paths.load_settings().get('auto_update', True) is not False


def current():
	"""What this copy is: {mode, branch, commit}."""
	if mode() == 'git':
		head, _ = _git('rev-parse', 'HEAD')
		br, _ = _git('rev-parse', '--abbrev-ref', 'HEAD')
		up, _ = _git('rev-parse', '--abbrev-ref', '--symbolic-full-name', '@{u}')
		return {'mode': 'git', 'branch': br, 'upstream': up, 'commit': head}
	rec = load_record()
	return {'mode': 'archive', 'branch': branch(), 'commit': rec.get('commit'), 'date': rec.get('date')}


# Copies without git ------------------------------------------------------------------------------------------------

def _archive_update(log=None):
	rec = load_record()
	br = branch()
	latest = _get(f'{API}/commits/{br}', accept='application/vnd.github.sha').decode().strip()
	if not re.fullmatch(r'[0-9a-f]{40}', latest):
		raise RuntimeError(f'unexpected answer of GitHub: {latest[:60]}')
	installed = rec.get('commit')
	if installed == latest:
		return {'updated': False, 'commit': latest, 'message': f'up to date ({br} {latest[:7]})'}
	if installed:
		cmp = json.loads(_get(f'{API}/compare/{installed}...{latest}'))
		if cmp.get('status') != 'ahead':
			return {'updated': False, 'commit': installed,
			        'message': f'the last commit of {br} ({latest[:7]}) is not newer than this copy ({installed[:7]}): '
			                   f'{cmp.get("status")}; left as it is'}
	_say(log, f'Updating eHabitat+ from GitHub ({REPO}, {br}: {(installed or "unknown")[:7]} -> {latest[:7]})...')
	data = _get(f'https://codeload.github.com/{REPO}/zip/{latest}', accept='application/zip', timeout=120)
	z = zipfile.ZipFile(io.BytesIO(data))
	members = {}
	for info in z.infolist():
		if info.is_dir() or '/' not in info.filename:
			continue
		rel = info.filename.split('/', 1)[1]
		if not rel or rel.startswith(NEVER) or '..' in rel.split('/'):
			continue
		members[rel] = info
	if REQUIRED not in members:
		raise RuntimeError(f'the commit {latest[:7]} of {br} has no app ({REQUIRED}): not applied')

	manifest = rec.get('files')
	new_manifest, updated, kept, removed = {}, [], [], []
	for rel, info in sorted(members.items()):
		content = z.read(info)
		h = _sha1(content)
		new_manifest[rel] = h
		dst = os.path.join(paths.ROOT, *rel.split('/'))
		if os.path.isfile(dst):
			cur = _file_sha1(dst)
			if cur == h:
				continue
			# changed by the user since it was installed (or, without a record, a file the user edits): kept
			changed = (cur != manifest[rel]) if manifest is not None and rel in manifest else \
				(manifest is not None or rel in PROTECTED)
			if changed:
				kept.append(rel)
				continue
		os.makedirs(os.path.dirname(dst), exist_ok=True)
		tmp = dst + '.ehabitat-tmp'
		with open(tmp, 'wb') as f:
			f.write(content)
		mode_bits = (info.external_attr >> 16) & 0o777
		if os.name != 'nt' and mode_bits:
			os.chmod(tmp, mode_bits)
		os.replace(tmp, dst)
		updated.append(rel)
	for rel, h in (manifest or {}).items():
		if rel in members:
			continue
		dst = os.path.join(paths.ROOT, *rel.split('/'))
		if os.path.isfile(dst) and _file_sha1(dst) == h:
			os.remove(dst)
			removed.append(rel)
	_save_record({'repo': REPO, 'branch': br, 'commit': latest,
	              'date': time.strftime('%Y-%m-%d %H:%M:%S'), 'files': new_manifest})
	env_changed = 'environment.yml' in updated
	if kept:
		_say(log, 'Kept with your changes: ' + ', '.join(kept))
	msg = f'updated to {br} {latest[:7]}: {len(updated)} files updated, {len(removed)} removed' if updated or removed \
		else f'up to date ({br} {latest[:7]})'
	_say(log, msg[0].upper() + msg[1:])
	return {'updated': bool(updated or removed), 'commit': latest, 'from': installed, 'files': updated,
	        'kept': kept, 'removed': removed, 'environment': env_changed, 'message': msg,
	        'app_changed': any(f.startswith('app/') for f in updated + removed)}


# Git clones --------------------------------------------------------------------------------------------------------

def _git_update(log=None):
	up, err = _git('rev-parse', '--abbrev-ref', '--symbolic-full-name', '@{u}')
	if not up:
		return {'updated': False, 'message': 'git clone without an upstream branch: not updated'}
	out, err = _git('fetch', '--quiet', timeout=60)
	if out is None:
		raise RuntimeError(f'git fetch failed: {err or "no answer"}')
	behind, _ = _git('rev-list', '--count', 'HEAD..@{u}')
	if not behind or behind == '0':
		head, _ = _git('rev-parse', 'HEAD')
		return {'updated': False, 'commit': head, 'message': f'up to date with {up}'}
	old, _ = _git('rev-parse', 'HEAD')
	_say(log, f'Updating eHabitat+ with git pull --ff-only ({up}, {behind} new commits)...')
	out, err = _git('merge', '--ff-only', '--quiet', '@{u}', timeout=120)
	if out is None:
		incoming, _ = _git('diff', '--name-only', 'HEAD', '@{u}')
		local, _ = _git('diff', '--name-only', 'HEAD')
		both = sorted(set((incoming or '').splitlines()) & set((local or '').splitlines()))
		why = f"local changes in {', '.join(both)}" if both else 'commits of its own'
		return {'updated': False, 'commit': old, 'message': f'{up} has {behind} new commits, but this git clone has '
		                                                     f'{why}: update it with git pull'}
	new, _ = _git('rev-parse', 'HEAD')
	files, _ = _git('diff', '--name-only', old, new)
	files = (files or '').splitlines()
	msg = f'updated to {up} {new[:7]} ({behind} commits)'
	_say(log, msg[0].upper() + msg[1:])
	return {'updated': True, 'commit': new, 'from': old, 'files': files, 'environment': 'environment.yml' in files,
	        'message': msg, 'app_changed': any(f.startswith('app/') for f in files)}


def update(log=None):
	"""Update this copy to the last version of GitHub. Returns {'updated', 'message', 'environment' (environment.yml
	changed), 'app_changed', 'files', ...}."""
	if any(j.kind in ('model', 'ecosim') or j.kind.startswith('gee') for j in jobs.running()):
		res = {'message': 'not updated while the model, the similarity or a download is running'}
	else:
		try:
			res = _git_update(log) if mode() == 'git' else _archive_update(log)
		except Exception as e:
			with _lock:
				_state.update(checked=time.time(), error=str(e))
			raise
	res = {'updated': False, 'environment': False, 'app_changed': False, 'files': [], **res}
	with _lock:
		_state.update(checked=time.time(), result=res, error=None)
		if res.get('updated'):
			_state['restart'] = True
	return res


def startup(log=None):
	"""Update before the window opens (quick: a few seconds). Returns the result, or None when skipped or offline."""
	if not enabled() or os.environ.get('EHAB_UPDATED'):
		return None
	try:
		res = update(log)
	except Exception as e:  # offline, GitHub unreachable, ...: the app opens as it is
		with _lock:
			_state.update(checked=time.time(), error=str(e))
		_say(log, f'Could not check for updates of eHabitat+ ({e}): opening this version')
		return None
	_say(log, f"eHabitat+: {res['message']}")
	return res


def restart_args():
	"""The command that opened this app (to open it again after an update)."""
	return [sys.executable, *(['-s'] if sys.flags.no_user_site else []), os.path.join(paths.APP_DIR, 'ehabitat.py'),
	        *sys.argv[1:]]


def reexec(res=None):
	"""Start this app again, with the new code, in this same process (before its window exists). RES: the update
	that was just made (its message is shown by the new process; with environment, it updates the environment)."""
	env = dict(os.environ, EHAB_UPDATED=(res or {}).get('message') or '1')
	if res and res.get('environment'):
		env['EHAB_ENV_UPDATE'] = '1'
	if os.name == 'nt':
		subprocess.Popen(restart_args(), env=env, cwd=paths.ROOT)
		os._exit(0)
	os.execve(sys.executable, restart_args(), env)


def spawn():
	"""Open a new instance of the app (the restart button: this one then closes its window)."""
	kw = {'creationflags': subprocess.DETACHED_PROCESS | NO_WINDOW} if os.name == 'nt' else {'start_new_session': True}
	subprocess.Popen(restart_args(), env=dict(os.environ, EHAB_UPDATED='1'), cwd=paths.ROOT,
	                 stdin=subprocess.DEVNULL, **kw)


# Environment -------------------------------------------------------------------------------------------------------

def _pip_requirements():
	"""The pip packages of environment.yml."""
	import yaml
	with open(os.path.join(paths.ROOT, 'environment.yml'), encoding='utf-8') as f:
		env = yaml.safe_load(f) or {}
	out = []
	for d in env.get('dependencies') or []:
		if isinstance(d, dict) and d.get('pip'):
			out += [str(x) for x in d['pip']]
	return out


WIN_CRAN = re.compile(r'\[win: CRAN ([A-Za-z0-9.]+)\]')  # mark of environment.yml: not in conda-forge for Windows
# the R packages of CRAN that R cannot load, installed as binary packages of Windows (as the installer of Windows does)
CRAN_INSTALL = ("p <- commandArgs(TRUE); m <- p[!sapply(p, requireNamespace, quietly = TRUE)]; if (length(m)) { "
                "message('From CRAN: ', paste(m, collapse = ', ')); install.packages(m, repos = 'https://cloud.r-project.org', "
                "type = 'binary'); m <- m[!sapply(m, requireNamespace, quietly = TRUE)]; if (length(m)) { "
                "message('R packages not installed: ', paste(m, collapse = ', ')); quit(status = 1) } }")


def conda_envfile():
	"""environment.yml as conda installs it here: on Windows, a temporary copy without the packages marked there with
	[win: CRAN name] (conda-forge does not have them for Windows). Returns (path, temporary, CRAN names)."""
	envfile = os.path.join(paths.ROOT, 'environment.yml')
	if os.name != 'nt':
		return envfile, False, []
	with open(envfile, encoding='utf-8') as f:
		lines = f.readlines()
	cran = [m.group(1) for m in map(WIN_CRAN.search, lines) if m]
	if not cran:
		return envfile, False, []
	fd, tmp = tempfile.mkstemp(prefix='ehabitat-environment-', suffix='.yml')
	with os.fdopen(fd, 'w', encoding='utf-8') as f:
		f.writelines(line for line in lines if not WIN_CRAN.search(line))
	return tmp, True, cran


def update_environment(job):
	"""Install what environment.yml asks for in this conda environment (micromamba or conda install -f), its pip
	packages that are missing and, on Windows, the R packages that come from CRAN there."""
	exes = [e for e in envsetup._mamba_candidates(sys.prefix) if os.path.basename(e).lower().startswith(('micromamba', 'conda', 'mamba'))]
	envfile, temporary, cran = conda_envfile()
	pkgs = None
	try:
		if exes:
			env = dict(os.environ)
			if os.path.basename(exes[0]).lower().startswith('micromamba'):
				env.setdefault('MAMBA_ROOT_PREFIX', os.path.dirname(exes[0]))
			if os.name == 'nt' and not env.get('CONDA_PKGS_DIRS'):
				# micromamba keeps the packages in pkgs\https\conda.anaconda.org\conda-forge\noarch\<package>: from the
				# folder of the app, the C++ headers of the compilers of R go over the 260 characters of the paths of
				# Windows, so the cache goes to a new short folder of the disk (as in the installer), deleted at the end.
				# mkdtemp makes it only for this user (Python >= 3.12.4 on Windows): a folder made in C:\ can be changed
				# by every user, and micromamba links the files of the cache into the environment
				try:
					pkgs = tempfile.mkdtemp(prefix='ehp-', dir=os.environ.get('SystemDrive', 'C:') + os.sep)
					env['CONDA_PKGS_DIRS'] = pkgs
				except OSError:
					pkgs = None
			code = jobs.run_command(job, [exes[0], 'install', '-y', '-p', sys.prefix, '-f', envfile], cwd=paths.ROOT, env=env)
			if code != 0:
				raise RuntimeError(f'the environment could not be updated (exit code {code})')
		else:
			job.say('micromamba or conda not found: only the pip packages are checked')
	finally:
		if temporary:
			os.remove(envfile)
		if pkgs:
			shutil.rmtree(pkgs, ignore_errors=True)
	rscript = shutil.which('Rscript') if cran else None
	if rscript:
		code = jobs.run_command(job, [rscript, '-e', CRAN_INSTALL, *cran], cwd=paths.ROOT)
		if code != 0:
			raise RuntimeError(f'the R packages of CRAN ({", ".join(cran)}) could not be installed')
	for req in _pip_requirements():
		name = re.split(r'[\s<>=!~@\[;]', req, 1)[0]
		try:
			importlib.metadata.version(name)
			continue
		except importlib.metadata.PackageNotFoundError:
			pass
		code = jobs.run_command(job, [sys.executable, '-m', 'pip', 'install', '--disable-pip-version-check', req],
		                        cwd=paths.ROOT)
		if code != 0:
			raise RuntimeError(f'{req} could not be installed')
	job.say('Environment updated: restart the app to use it.')
	with _lock:
		_state['restart'] = True
	return True


def start_environment_update():
	job = jobs.start_function('update-env', 'Update of the environment', update_environment)
	with _lock:
		_state['env_job'] = job.id
	return job.id


def status():
	with _lock:
		st = dict(_state)
	st.update(current())
	st['auto'] = enabled()
	st['repo'] = REPO
	just = os.environ.get('EHAB_UPDATED')
	st['just_updated'] = just if just and just != '1' else None  # the update made when this app opened
	if not st['result'] and st['just_updated']:
		st['result'] = {'updated': True, 'message': st['just_updated']}
	return st
