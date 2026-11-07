#!/usr/bin/env python3
#### eHabitat+ app: window with the inputs, the HFT model, MonitorEO and the results
#### License: CC BY-SA 3.0

"""
Open the app with the Python of the environment of eHabitat+ (environment.yml):

    python app/ehabitat.py            # in its window (pywebview)
    python app/ehabitat.py --browser  # in the web browser, without pywebview

The installers of app/installer create that environment and a shortcut on the Desktop (see the page
Desktop app of the wiki: https://github.com/javimarlop/eHabpy/wiki/Desktop-app).
The app works on the copy of eHabpy it belongs to: inVars, segm_optim (config.py, pa/, results/), ...
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if sys.stdout is None:  # pythonw.exe (the shortcut of Windows), without a console: the messages go to ehabitat.log,
	# next to the environment made by the installer (<folder of eHabitat+>\env), or in ~/.ehabitat
	_dirs = [os.path.dirname(sys.prefix)] if os.path.basename(sys.prefix) == 'env' else []
	for _logdir in _dirs + [os.path.join(os.path.expanduser('~'), '.ehabitat')]:
		try:
			os.makedirs(_logdir, exist_ok=True)
			sys.stdout = sys.stderr = open(os.path.join(_logdir, 'ehabitat.log'), 'a', encoding='utf-8',
			                               errors='replace', buffering=1)
			break
		except OSError:  # a folder this user cannot write: the next one or, without any, no messages (as before)
			pass

from ehabitat_app.main import main  # noqa: E402

if __name__ == '__main__':
	main()
