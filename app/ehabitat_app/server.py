"""Local HTTP server (127.0.0.1, a random port and a secret token in every path) of the interface (ui/), of the files
of the copy of eHabpy that the interface shows (plots of the results) and, when the app runs in a web browser instead
of its window, of the functions of the API (POST <token>/api/<function>, with the arguments as JSON)."""

import json
import mimetypes
import os
import secrets
import threading
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from . import paths

FILE_TYPES = ('.png', '.jpg', '.jpeg', '.svg', '.csv', '.txt', '.json', '.geojson', '.pdf', '.yaml')


class Server:
	def __init__(self, api):
		self.api = api
		self.token = secrets.token_hex(16)
		server = self

		class Handler(BaseHTTPRequestHandler):
			def log_message(self, fmt, *args):
				if os.environ.get('EHAB_HTTP_LOG'):  # each request, to debug the interface
					print('http: ' + (fmt % args).replace(server.token, '<token>'), flush=True)

			def _host_ok(self):
				# against DNS rebinding: only requests to this address
				return self.headers.get('Host') in (f'127.0.0.1:{server.port}', f'localhost:{server.port}')

			def _parts(self):
				path = urllib.parse.urlparse(self.path).path
				parts = path.lstrip('/').split('/', 2)
				if len(parts) < 2 or not secrets.compare_digest(parts[0], server.token):
					return None
				return parts[1], urllib.parse.unquote(parts[2]) if len(parts) > 2 else ''

			def _send(self, code, body, ctype='application/json; charset=utf-8', cache=False):
				self.send_response(code)
				self.send_header('Content-Type', ctype)
				self.send_header('Content-Length', str(len(body)))
				self.send_header('Cache-Control', 'max-age=3600' if cache else 'no-store')
				self.send_header('X-Content-Type-Options', 'nosniff')
				self.end_headers()
				self.wfile.write(body)

			def do_GET(self):
				parts = self._parts() if self._host_ok() else None
				if not parts:
					return self._send(404, b'{}')
				area, rest = parts
				if area == 'ui':
					base = paths.UI
					allowed = None
				elif area == 'files':
					base = paths.ROOT
					allowed = FILE_TYPES
				else:
					return self._send(404, b'{}')
				target = os.path.normpath(os.path.join(base, rest or 'index.html'))
				if not paths.inside(target, base) or not os.path.isfile(target):
					return self._send(404, b'{}')
				if allowed and not target.lower().endswith(allowed):
					return self._send(403, b'{}')
				ctype = mimetypes.guess_type(target)[0] or 'application/octet-stream'
				if ctype.startswith('text/') or ctype in ('application/javascript', 'application/json'):
					ctype += '; charset=utf-8'
				with open(target, 'rb') as f:
					body = f.read()
				self._send(200, body, ctype, cache=(area == 'files' and '?v=' in self.path))

			def do_POST(self):
				parts = self._parts() if self._host_ok() else None
				if not parts or parts[0] != 'api':
					return self._send(404, b'{}')
				name = parts[1]
				fn = getattr(server.api, name, None)
				if name.startswith('_') or not callable(fn):
					return self._send(404, json.dumps({'ok': False, 'error': f'no function {name}'}).encode())
				n = int(self.headers.get('Content-Length') or 0)
				try:
					args = json.loads(self.rfile.read(n) or b'[]') if n else []
				except ValueError:
					return self._send(400, b'{}')
				res = fn(*args) if isinstance(args, list) else fn(**args)
				self._send(200, json.dumps(res, default=str).encode())

		self.httpd = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
		self.httpd.daemon_threads = True
		self.port = self.httpd.server_address[1]

	def start(self):
		threading.Thread(target=self.httpd.serve_forever, daemon=True, name='http').start()
		return self

	def url(self, path=''):
		return f'http://127.0.0.1:{self.port}/{self.token}/{path}'

	def stop(self):
		self.httpd.shutdown()
