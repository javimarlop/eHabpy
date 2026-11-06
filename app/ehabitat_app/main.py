"""Open the app: activate the environment, start the local server and show the interface in its window (pywebview) or,
with --browser or without pywebview, in the web browser."""

import argparse
import os
import sys
import time
import webbrowser

from . import __version__, envsetup, jobs


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
	if not args.browser:
		try:
			import webview
		except ImportError as e:
			print(f'pywebview is not installed ({e}): opening the interface in the web browser', flush=True)
	if webview is None:
		print(f'Interface: {url}\nPress Ctrl+C to close the app.', flush=True)
		if not args.no_open:
			webbrowser.open(url)
		try:
			while True:
				time.sleep(3600)
		except KeyboardInterrupt:
			pass
		finally:
			jobs.cancel_all()
		return

	window = webview.create_window('eHabitat+', url, js_api=api, width=1440, height=920, min_size=(1024, 680),
	                               text_select=True)
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
	finally:
		jobs.cancel_all()
