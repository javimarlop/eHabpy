"""Open the app: activate the environment, start the local server and show the interface in its window (pywebview) or,
with --browser or without pywebview, in the web browser."""

import argparse
import os
import subprocess
import sys
import time
import webbrowser

from . import __version__, envsetup, jobs, paths


def main():
	parser = argparse.ArgumentParser(description='eHabitat+: inputs, HFT model, MonitorEO and results in a window.')
	parser.add_argument('--browser', action='store_true', help='open the interface in the web browser, not in a window')
	parser.add_argument('--debug', action='store_true', help='with the developer tools of the window')
	parser.add_argument('--no-open', action='store_true', help='with --browser, only print the address (tests)')
	args = parser.parse_args()
	try:
		sys.stdout.reconfigure(line_buffering=True)
	except (AttributeError, ValueError):
		pass

	activation = envsetup.activate()
	# The last version of GitHub before anything else is imported; with new code of the app, start again with it
	from . import updater
	res = updater.startup()
	if res and res.get('updated') and res.get('app_changed'):
		print('Opening the new version...', flush=True)
		updater.reexec(res)
	from .api import Api
	from .server import Server
	api = Api()
	api._activation = activation
	server = Server(api).start()
	url = server.url('ui/index.html')
	print(f'eHabitat+ app {__version__} | environment: {sys.prefix} ({activation})', flush=True)
	if (res and res.get('environment')) or os.environ.get('EHAB_ENV_UPDATE'):
		os.environ.pop('EHAB_ENV_UPDATE', None)
		updater.start_environment_update()  # environment.yml changed: in the background, with the window open

	webview = None
	why = None if args.browser else linux_window_problem()
	if why:
		print(f'The window cannot be opened ({why}): opening the interface in the web browser', flush=True)
	elif not args.browser:
		try:
			import webview
		except ImportError as e:
			print(f'pywebview is not installed ({e}): opening the interface in the web browser', flush=True)
	if webview is None:
		return in_browser(url, args.no_open)

	try:
		window = webview.create_window('eHabitat+', url, width=1440, height=920, min_size=(1024, 680), text_select=True)
	except Exception as e:
		print(f'The window cannot be created ({e}): opening the interface in the web browser', flush=True)
		return in_browser(url, args.no_open)
	api._window = window

	def closing():
		running = jobs.running()
		if not running:
			return True
		titles = ', '.join(j.title for j in running)
		try:
			ok = window.create_confirmation_dialog('Close eHabitat+', f'There are jobs running ({titles}). '
			                                       'Closing the app stops them. Close it?')
		except Exception:
			ok = True
		return bool(ok)

	try:
		window.events.closing += closing
	except Exception:
		pass
	try:
		webview.start(debug=args.debug, private_mode=False,
		              storage_path=os.path.join(os.path.expanduser('~'), '.ehabitat', 'webview'))
	except Exception as e:
		# no GUI toolkit for the window (e.g. Qt or GTK missing on Linux): the same interface in the web browser
		print(f'The window cannot be opened ({type(e).__name__}: {e}): opening the interface in the web browser',
		      flush=True)
		api._window = None
		return in_browser(url, args.no_open)
	finally:
		if api._window is not None:
			jobs.cancel_all()


def linux_window_problem():
	"""On Linux, why the window (Qt WebEngine) cannot open, or None: Qt aborts the whole process when it cannot reach
	the display or load its platform plugin, so it is tried first in another process."""
	if not sys.platform.startswith('linux'):
		return None
	if not (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY')):
		return 'no graphical display'
	code = ('import sys; from qtpy import QtWebEngineWidgets; from qtpy.QtWidgets import QApplication; '
	        'QApplication(sys.argv[:1])')
	try:
		r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=90,
		                   stdin=subprocess.DEVNULL)
	except (OSError, subprocess.SubprocessError) as e:
		return str(e)
	if r.returncode != 0:
		lines = [x for x in (r.stderr or '').strip().splitlines() if x.strip()]
		return 'Qt: ' + (lines[-1] if lines else f'exit code {r.returncode}')
	return None


def in_browser(url, no_open=False):
	"""Show the interface in the web browser and serve it until Ctrl+C."""
	print(f'Interface: {url}\nPress Ctrl+C to close the app.', flush=True)
	if not no_open and not (paths.wsl() and paths.open_in_windows(url)):  # WSL: the web browser of Windows
		webbrowser.open(url)
	try:
		while True:
			time.sleep(3600)
	except KeyboardInterrupt:
		pass
	finally:
		jobs.cancel_all()
