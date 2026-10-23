#### eHabitat+ HFT module: runs the steps that are out of date, in order
#### License: CC BY-SA 3.0

"""
Run from the segm_optim folder, after python3 setup.py:

    python3 run_hft.py               # run every step that is needed
    python3 run_hft.py --status      # only show what would run and why
    python3 run_hft.py --from stats  # run again from a step (segm, stats, moranvar or plots)

A step is needed when it was never run, did not finish, or when something it
depends on changed since it ran: its options in config.py, palist.csv, the
input files or its code (see hft_state.py). Then the steps after it are run
too. When the input files changed since setup.py imported them, setup.py is
run first with the options it was run with.
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


def setup_command(rec):
	cmd = [sys.executable, 'setup.py', '--grassdb', rec['grassdb'], '--col-id', rec['col_id']]
	if rec.get('study_area_file'):
		cmd += ['--study-area', rec['study_area_file']]
	if rec.get('layer'):
		cmd += ['--layer', rec['layer']]
	return cmd


def main():
	parser = argparse.ArgumentParser(description='Run the steps of the HFT module that are out of date.')
	parser.add_argument('--status', action='store_true', help='only show what would run and why')
	parser.add_argument('--from', dest='first', choices=hft_state.NAMES,
	                    help='run again from this step, even if it is up to date')
	args = parser.parse_args()

	cfg, state = hft_state.read_config(), hft_state.load()
	project = os.path.join(os.path.expanduser(str(cfg.get('GRASSDB'))), str(cfg.get('GRASSLOC')), 'PERMANENT')
	if not os.path.isdir(project):
		sys.exit(f'ERROR: GRASS project not found: {project}\nRun python3 setup.py first (see README.md)')

	hft_state.print_status(cfg, state)
	changes = hft_state.setup_changes(cfg, state)
	todo = [name for name, r in hft_state.plan(cfg, state) if r]
	if changes:
		todo = list(hft_state.NAMES)  # importing the data again makes every step out of date
	if args.first:
		todo = hft_state.NAMES[hft_state.NAMES.index(args.first):]
	if not todo and not changes:
		print(f"\nEverything is up to date: the results are in results/{cfg.get('STUDY_AREA')}")
		return
	# moranvar.py runs moranvar_plots.R at the end
	commands = [hft_state.COMMAND[n] for n in todo if not (n == 'plots' and 'moranvar' in todo)]
	if changes:
		commands.insert(0, 'python3 setup.py (imports the files that changed)')
	print('\nTo run: ' + ', '.join(commands))
	if args.status:
		return

	env = dict(os.environ, PYTHONNOUSERSITE='1', EHAB_PYTHON=sys.executable)
	try:  # as 'ulimit -n 8192'
		import resource
		soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
		want = 8192 if hard == resource.RLIM_INFINITY else min(8192, hard)
		if soft != resource.RLIM_INFINITY and soft < want:
			resource.setrlimit(resource.RLIMIT_NOFILE, (want, hard))
	except (ImportError, ValueError, OSError):
		pass

	if changes:
		run(setup_command(hft_state.get_setup(cfg.get('STUDY_AREA'), state)), env)
	if 'segm' in todo and args.first == 'segm':
		env['EHAB_FORCE_RESTART'] = '1'  # same inputs: it would otherwise find nothing to do
	for name in todo:
		if name == 'plots' and 'moranvar' in todo:
			continue
		cmd = hft_state.COMMAND[name].split()
		run([sys.executable] + cmd[1:] if cmd[0] == 'python3' else cmd, env)

	print(f"\nDone: the results are in results/{cfg.get('STUDY_AREA')}")
	hft_state.print_status()


if __name__ == '__main__':
	os.chdir(hft_state.HERE)
	main()
