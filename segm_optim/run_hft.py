#### eHabitat+ HFT module: runs the steps that are out of date, in order
#### License: CC BY-SA 3.0

"""
Run from the segm_optim folder, after python3 setup.py:

    python3 run_hft.py                     # run every step that is needed
    python3 run_hft.py --status            # only show what would run and why
    python3 run_hft.py --scenario NAME     # results in results/<study area>/NAME, without asking
    python3 run_hft.py --from stats        # run again from a step (segm, stats, moranvar or plots)

The results are written to results/<STUDY_AREA>/<SCENARIO> (SCENARIO in
config.py or, if empty, the name of INVARS_FOLDER). When run in a terminal, it
shows the scenario and asks whether to change it before starting, and it warns
when the scenario already has results made with other options (they would be
moved to <SCENARIO>_previous).

A step is needed when it was never run, did not finish, or when something it
depends on changed since it ran: its options in config.py, palist.csv, the
input files or its code (see hft_state.py). Then the steps after it are run
too. When the input files changed since setup.py imported them, or
INVARS_FOLDER changed, setup.py is run first with the options it was run with.
"""

import argparse
import os
import subprocess
import sys

import hft_state


def run(cmd, env):
	print(f"\n>>> {' '.join(cmd)}", flush=True)
	r = subprocess.run(cmd, cwd=hft_state.HERE, env=env)
	if r.returncode != 0:
		sys.exit(f"ERROR: {' '.join(cmd)} failed (exit code {r.returncode}): see the messages above")


def setup_command(rec, cfg):
	"""setup.py with the options of its last run (those of config.py if it left no record)."""
	rec = rec or {}
	cmd = [sys.executable, 'setup.py', '--grassdb', rec.get('grassdb', str(cfg.get('GRASSDB'))),
	       '--col-id', rec.get('col_id', str(cfg.get('COL_ID')))]
	if rec.get('study_area_file'):
		cmd += ['--study-area', rec['study_area_file']]
	if rec.get('layer'):
		cmd += ['--layer', rec['layer']]
	return cmd


def ask_name(prompt, default):
	"""Scenario name typed by the user (DEFAULT if only Enter is pressed)."""
	while True:
		try:
			name = input(prompt).strip() or default
		except EOFError:
			return default
		err = hft_state.scenario_error(name)
		if not err:
			return name
		print(f"'{name}' is not a valid scenario name: {err}")


def choose_scenario(args, cfg, state, interactive, compare=True):
	"""Scenario of this run: --scenario NAME, or the one of config.py, which can be changed when
	run in a terminal. Warns when the scenario already has results made with other options (unless
	not COMPARE: setup.py has still to fill config.py with the variables to use)."""
	name = hft_state.scenario_name(cfg)
	if args.scenario is not None:
		name = args.scenario
		err = hft_state.scenario_error(name)
		if err:
			sys.exit(f"ERROR: '{name}' is not a valid scenario name: {err}")
	elif interactive:
		print(f"Scenario: {name} (results in {hft_state.results_dir(cfg)})")
		name = ask_name('Press Enter to keep it, or type another scenario name: ', name)
	err = hft_state.scenario_error(name)
	if err:
		sys.exit(f"ERROR: '{name}' is not a valid scenario name (SCENARIO or INVARS_FOLDER in config.py): {err}")
	while compare:
		this = dict(cfg, SCENARIO=name)
		folder = hft_state.results_dir(this)
		changes = hft_state.scenario_changes(this, state)
		if changes == [] or not (os.path.isdir(folder) and os.listdir(folder)):
			return name
		made = f"with other options ({'; '.join(changes)})" if changes else 'with options that are not recorded'
		print(f"\nWARNING: {folder} already has results made {made}.\n"
		      f"Running again moves them to {folder}_previous (replacing what is there).")
		if not interactive:
			return name
		new = ask_name('Press Enter to replace them, or type another scenario name: ', name)
		if new == name:
			return name
		name = new
	return name


def main():
	parser = argparse.ArgumentParser(description='Run the steps of the HFT module that are out of date.')
	parser.add_argument('--status', action='store_true', help='only show what would run and why')
	parser.add_argument('--scenario', metavar='NAME',
	                    help='scenario of the results (results/<study area>/NAME), without asking; '
	                         'it is saved as SCENARIO in config.py')
	parser.add_argument('--from', dest='first', choices=hft_state.NAMES,
	                    help='run again from this step, even if it is up to date')
	args = parser.parse_args()

	cfg, state = hft_state.read_config(), hft_state.load()
	project = os.path.join(os.path.expanduser(str(cfg.get('GRASSDB'))), str(cfg.get('GRASSLOC')), 'PERMANENT')
	if not os.path.isdir(project):
		sys.exit(f'ERROR: GRASS project not found: {project}\nRun python3 setup.py first (see README.md)')

	env = dict(os.environ, PYTHONNOUSERSITE='1', EHAB_PYTHON=sys.executable)
	try:  # as 'ulimit -n 8192'
		import resource
		soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
		want = 8192 if hard == resource.RLIM_INFINITY else min(8192, hard)
		if soft != resource.RLIM_INFINITY and soft < want:
			resource.setrlimit(resource.RLIMIT_NOFILE, (want, hard))
	except (ImportError, ValueError, OSError):
		pass

	# setup.py first when the input data or INVARS_FOLDER changed: it fills config.py with the variables
	# to use, which the scenario is compared with (with --status it is only shown)
	changes = hft_state.setup_changes(cfg, state)
	if changes and not args.status:
		print(f"The input data changed: {'; '.join(changes)}")
		run(setup_command(hft_state.get_setup(cfg.get('STUDY_AREA'), state), cfg), env)
		cfg, state = hft_state.read_config(), hft_state.load()

	# The scenario is asked only in a terminal (and not with --status)
	scenario = choose_scenario(args, cfg, state, interactive=sys.stdin.isatty() and not args.status,
	                           compare=not (changes and args.status))
	if args.status:
		cfg = dict(cfg, SCENARIO=scenario)  # shown, not saved
	elif scenario != hft_state.scenario_name(cfg):
		hft_state.write_option('SCENARIO', repr(scenario))
		print(f'config.py: SCENARIO = {scenario!r}')
		cfg = hft_state.read_config()
	print()

	hft_state.print_status(cfg, state)
	todo = [name for name, r in hft_state.plan(cfg, state) if r]
	if changes:
		todo = list(hft_state.NAMES)  # importing the data again makes every step out of date
	if args.first:
		todo = hft_state.NAMES[hft_state.NAMES.index(args.first):]
	if not todo and not changes:
		print(f'\nEverything is up to date: the results are in {hft_state.results_dir(cfg)}')
		return
	# moranvar.py runs moranvar_plots.R at the end
	commands = [hft_state.COMMAND[n] for n in todo if not (n == 'plots' and 'moranvar' in todo)]
	if changes and args.status:
		commands.insert(0, 'python3 setup.py (imports the variables to use)')
	print('\nTo run: ' + ', '.join(commands))
	if args.status:
		return

	if 'segm' in todo and args.first == 'segm':
		env['EHAB_FORCE_RESTART'] = '1'  # same inputs: it would otherwise find nothing to do
	for name in todo:
		if name == 'plots' and 'moranvar' in todo:
			continue
		cmd = hft_state.COMMAND[name].split()
		run([sys.executable] + cmd[1:] if cmd[0] == 'python3' else cmd, env)

	print(f'\nDone: the results are in {hft_state.results_dir(cfg)}')
	hft_state.print_status()


if __name__ == '__main__':
	os.chdir(hft_state.HERE)
	main()
