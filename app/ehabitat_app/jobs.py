"""Background jobs: the scripts of eHabpy (one after the other) or a Python function, with their log (read by the
interface while they run), their progress and their cancellation (the whole tree of processes)."""

import itertools
import os
import signal
import subprocess
import sys
import threading
import time
import traceback

from . import paths

MAX_LINES = 20000
NO_WINDOW = getattr(subprocess, 'CREATE_NO_WINDOW', 0)
NEW_GROUP = getattr(subprocess, 'CREATE_NEW_PROCESS_GROUP', 0)


class Cancelled(Exception):
	pass


class Job:
	def __init__(self, kind, title):
		self.id = next(_ids)
		self.kind, self.title = kind, title
		self.state = 'running'  # running, done, failed, cancelled
		self.lines = []
		self.dropped = 0  # first lines removed when the log is longer than MAX_LINES
		self.overwrite = False  # the next text replaces the last line (\r of the progress bars)
		self.backspaces = 0  # the next text replaces these last characters of the last line
		self.progress = {}
		self.result = None
		self.error = None
		self.started, self.ended = time.time(), None
		self.proc = None
		self.cancel_requested = False
		self._lock = threading.Lock()

	# log ------------------------------------------------------------------------------------------------------------
	def log(self, text):
		"""Add TEXT to the log: each line; a carriage return makes the next text replace the last line, and N backspaces
		(the percentages of GRASS) make it replace the last N characters."""
		with self._lock:
			for part in _split(str(text)):
				if part == '\n':
					self.overwrite = False
					self.backspaces = 0
					if not self.lines:
						self.lines.append('')
					self.lines.append('')
				elif part == '\r':
					self.overwrite = True
					self.backspaces = 0
				elif part == '\b':
					self.backspaces += 1
				else:
					if self.backspaces and self.lines:
						self.lines[-1] = self.lines[-1][:-self.backspaces] if self.backspaces < len(self.lines[-1]) else ''
						self.backspaces = 0
					if self.overwrite and self.lines:
						self.lines[-1] = part
						self.overwrite = False
					elif self.lines:
						self.lines[-1] += part
					else:
						self.lines.append(part)
			if len(self.lines) > MAX_LINES:
				n = len(self.lines) - MAX_LINES
				del self.lines[:n]
				self.dropped += n

	def say(self, text):
		"""Add TEXT as whole lines."""
		self.log(str(text).rstrip('\n') + '\n')

	def check(self):
		if self.cancel_requested:
			raise Cancelled()

	def snapshot(self, since=0):
		"""State and the lines from SINCE (number of the line, counting the dropped ones), resending the last line the
		reader has, which may have changed: the reader replaces its lines from 'start'."""
		with self._lock:
			total = self.dropped + len(self.lines)
			start = max(self.dropped, min(int(since or 0), total) - 1)
			return {'id': self.id, 'kind': self.kind, 'title': self.title, 'state': self.state,
			        'start': start, 'lines': self.lines[start - self.dropped:], 'total': total,
			        'progress': self.progress, 'result': self.result if self.state != 'running' else None,
			        'error': self.error, 'elapsed': round((self.ended or time.time()) - self.started, 1)}

	def finish(self, state, error=None):
		with self._lock:
			if self.state == 'running':
				self.state = state
				self.error = error
				self.ended = time.time()


def _split(text):
	out, buf = [], ''
	for ch in text:
		if ch in '\r\n\b':
			if buf:
				out.append(buf)
				buf = ''
			out.append(ch)
		else:
			buf += ch
	if buf:
		out.append(buf)
	# \r\n is a new line
	merged = []
	for p in out:
		if p == '\n' and merged and merged[-1] == '\r':
			merged[-1] = '\n'
		else:
			merged.append(p)
	return merged


_ids = itertools.count(1)
_jobs = {}
_jobs_lock = threading.Lock()


def get(job_id):
	return _jobs.get(int(job_id))


def running(kind=None):
	return [j for j in _jobs.values() if j.state == 'running' and (kind is None or j.kind == kind)]


def _register(job):
	with _jobs_lock:
		_jobs[job.id] = job
		# keep the last 30 finished jobs
		old = sorted((j for j in _jobs.values() if j.state != 'running'), key=lambda j: j.id)[:-30]
		for j in old:
			_jobs.pop(j.id, None)


def start_function(kind, title, fn, *args, **kwargs):
	"""Run FN(job, *args, **kwargs) in a thread; what it returns is the result of the job."""
	job = Job(kind, title)
	_register(job)

	def work():
		try:
			job.result = fn(job, *args, **kwargs)
			job.finish('cancelled' if job.cancel_requested else 'done')
		except Cancelled:
			job.say('Cancelled.')
			job.finish('cancelled')
		except SystemExit as e:  # fail() of the scripts of eHabpy
			msg = str(e.code) if e.code not in (None, 0) else 'stopped'
			job.say(msg)
			job.finish('failed', msg)
		except Exception as e:
			job.say(''.join(traceback.format_exception_only(type(e), e)).strip())
			job.finish('failed', str(e))

	threading.Thread(target=work, daemon=True, name=f'job-{job.id}').start()
	return job


def run_command(job, cmd, cwd=paths.SEGM, env=None):
	"""Run CMD writing its output to the log of JOB; returns its exit code. Cancelling the job stops it (and the
	processes it started)."""
	job.check()
	shown = [os.path.basename(c) if i == 0 else c for i, c in enumerate(cmd)]
	job.say(f"\n>>> {' '.join(shown)}")
	env = dict(os.environ if env is None else env)
	env.setdefault('PYTHONUNBUFFERED', '1')
	env.setdefault('PYTHONIOENCODING', 'utf-8')
	kwargs = {}
	if os.name == 'nt':
		kwargs['creationflags'] = NO_WINDOW | NEW_GROUP
	else:
		kwargs['start_new_session'] = True
	proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
	                        stderr=subprocess.STDOUT, **kwargs)
	job.proc = proc
	try:
		fd = proc.stdout.fileno()
		decoder = _decoder()
		while True:
			chunk = os.read(fd, 8192)
			if not chunk:
				break
			job.log(decoder.decode(chunk))
		job.log(decoder.decode(b'', final=True))
		code = proc.wait()
	finally:
		job.proc = None
		proc.stdout.close()
	job.log('\n' if job.lines and job.lines[-1] else '')
	if job.cancel_requested:
		raise Cancelled()
	return code


def _decoder():
	import codecs
	return codecs.getincrementaldecoder('utf-8')(errors='replace')


def start_commands(kind, title, commands, env=None, cwd=paths.SEGM, after=None, intro=None):
	"""Run COMMANDS one after the other, stopping at the first that fails; AFTER(job) runs at the end if all ran.
	INTRO: text written first to the log."""
	def work(job):
		if intro:
			job.say(intro)
		for cmd in commands:
			code = run_command(job, cmd, cwd=cwd, env=env)
			if code != 0:
				raise RuntimeError(f"{os.path.basename(cmd[1] if len(cmd) > 1 and cmd[0] == sys.executable else cmd[0])} "
				                   f"failed (exit code {code}): see the messages above")
		return after(job) if after else None
	return start_function(kind, title, work)


def cancel(job_id):
	job = get(job_id)
	if not job or job.state != 'running':
		return False
	job.cancel_requested = True
	proc = job.proc
	if proc and proc.poll() is None:
		job.say('\nStopping...')
		kill_tree(proc)
	return True


def kill_tree(proc):
	"""Stop PROC and every process it started."""
	try:
		if os.name == 'nt':
			subprocess.run(['taskkill', '/T', '/F', '/PID', str(proc.pid)], capture_output=True, creationflags=NO_WINDOW)
		else:
			os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
			for _ in range(50):
				if proc.poll() is not None:
					return
				time.sleep(0.1)
			os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
	except (OSError, ProcessLookupError):
		pass


def cancel_all():
	for j in running():
		cancel(j.id)
