"""The changelog of eHabitat+: the section "## Changelog" at the end of README.md, what each release brings, newest
first (as that of MonitorEO, https://github.com/monitoreo-obsnev/monitoreo#changelog).

    python .github/changelog.py check                 # its format, and that it has the releases of the tags
    python .github/changelog.py notes 1.5             # the entries of a release in Markdown (its GitHub release notes)
    python .github/changelog.py unrecorded OLD NEW    # whether the commits OLD..NEW change eHabitat+ without an entry
    python .github/changelog.py release [X.Y.Z]       # makes the release of "### Unreleased" (commit and tag)
    python .github/changelog.py pending               # whether "### Unreleased" has entries (pushed without a release)

Each release is a heading with its tag and the date of its GitHub release, followed by its entries, one per line, each
with the parts of eHabitat+ it reaches (AREAS) and "(fix)" for the fixes:

    ### [1.4](https://github.com/javimarlop/eHabpy/releases/tag/1.4) · 2026-10-27

    - **HFT model**: New CLASS_METHOD option in config.py ...
    - **Similarity** (fix): The similarity indices are right when a mask is used.

Above the releases, "### Unreleased" holds the entries of the changes not yet released (the desktop app has those of
master already: it updates itself from master every time it opens). Every change that users notice adds its entry
there, in the same commit, and every push to master with entries is a release: "release" renames "### Unreleased" to
the next version (the last number plus one: 1.6.1 after 1.6, then 1.6.2...; or the one given, e.g. 1.7) with the date
of today in UTC (that of its GitHub release), commits README.md alone ("Release X.Y.Z") and tags that commit (an
annotated tag). Pushed with

    git push --atomic --follow-tags origin master

the tag makes .github/workflows/release.yml create its GitHub release (notes: "notes X.Y.Z"), which Zenodo archives.
If the push is rejected (master moved), undo the release (git tag -d X.Y.Z; git reset --keep HEAD~1), pull and make
it again: a heading of a release without its tag is a problem for "check".
"--trailer TEXT" (as often as needed) adds a trailer to the commit of the release. Only the standard library: it runs
in GitHub Actions (.github/workflows/changelog.yml and release.yml).
"""

import datetime
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = 'https://github.com/javimarlop/eHabpy'
SECTION = '## Changelog'
UNRELEASED = '### Unreleased'
AREAS = ('eHabitat+', 'Desktop app', 'Installers', 'HFT model', 'Similarity', 'Earth Engine variables', 'Documentation')
RELEASE = re.compile(r'^### \[(?P<tag>\d+(?:\.\d+)+)\]\(' + re.escape(REPO) + r'/releases/tag/(?P<link>[^)]+)\)'
                     r' · (?P<date>\d{4}-\d{2}-\d{2})$')
ENTRY = re.compile(r'^- \*\*(?P<areas>[^*]+)\*\*(?P<fix> \(fix\))?: (?P<text>\S.*)$')
FORMAT = (f'a release is "### [<tag>]({REPO}/releases/tag/<tag>) · <YYYY-MM-DD>", and an entry is '
          '"- **<areas>**: <text>." or "- **<areas>** (fix): <text>."')
# The files of eHabitat+ that users notice when they change (code, environment, installers); the others (README, wiki
# images, tests of the maintainers...) need no entry
CODE = re.compile(r'^(app/|segm_optim/.*\.(py|R|yaml)$|[^/]+\.(py|R)$|environment\.yml$|requirements\.txt$|'
                  r'conf_grass7eHabplus\.sh$)')


def version(v):
	"""A version as a tuple to compare: 1.4 = (1, 4, 0)."""
	parts = [int(p) for p in str(v).strip().lstrip('vV').split('.')]
	return tuple(parts + [0] * (3 - len(parts)))


def section(text):
	"""The lines of the section "## Changelog" (without its heading), with their line numbers in README.md."""
	lines = text.splitlines()
	starts = [i for i, line in enumerate(lines) if line.strip() == SECTION]
	if len(starts) != 1:
		raise ValueError(f'README.md has {len(starts)} sections "{SECTION}" (one expected)')
	out = []
	for i in range(starts[0] + 1, len(lines)):
		if lines[i].startswith('## ') or lines[i].strip() == '---':
			break
		out.append((i + 1, lines[i]))
	return out


def parse(text):
	"""The changelog of README.md: {'unreleased': entries or None, 'releases': [{tag, date, entries, line}]}, newest
	first. An entry is {areas, fix, text, line}. Raises ValueError, with the line, if something is malformed."""
	unreleased, releases, current = None, [], None
	for n, line in section(text):
		if line.startswith('### '):
			if line.strip() == UNRELEASED:
				if releases or unreleased is not None:
					raise ValueError(f'README.md:{n}: "{UNRELEASED}" goes once, above the releases')
				unreleased = current = []
				continue
			m = RELEASE.match(line.rstrip())
			if not m:
				raise ValueError(f'README.md:{n}: a heading of the changelog that is not a release: {FORMAT}')
			if m['link'] != m['tag']:
				raise ValueError(f'README.md:{n}: the link of {m["tag"]} goes to the release {m["link"]}')
			current = []
			releases.append({'tag': m['tag'], 'date': m['date'], 'entries': current, 'line': n})
		elif line.startswith('- ') or (current is not None and line.strip()):
			if current is None:
				raise ValueError(f'README.md:{n}: an entry above the first release')
			m = ENTRY.match(line.rstrip())
			if not m:
				raise ValueError(f'README.md:{n}: not an entry of the changelog (one per line): {FORMAT}')
			areas = [a.strip() for a in m['areas'].split(',')]
			wrong = [a for a in areas if a not in AREAS]
			if wrong:
				raise ValueError(f'README.md:{n}: unknown area {", ".join(wrong)} (the areas: {", ".join(AREAS)})')
			current.append({'areas': areas, 'fix': bool(m['fix']), 'text': m['text'].strip(), 'line': n})
	return {'unreleased': unreleased, 'releases': releases}


def git(*args):
	"""The output of git in the repository, or None if it cannot run."""
	try:
		r = subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True, timeout=60)
	except (OSError, subprocess.SubprocessError):
		return None
	return r.stdout if r.returncode == 0 else None


def tags():
	"""The release tags of the repository (1.0, 1.4, ...), or None without git."""
	out = git('tag', '--list')
	if out is None:
		return None
	return [t for t in out.split() if re.fullmatch(r'v?\d+(?:\.\d+)+', t)]


def problems(text, known_tags=None):
	"""What is wrong with the changelog of README.md, against the release tags (by default, those of git): a list of
	texts, empty if nothing."""
	try:
		log = parse(text)
	except ValueError as e:
		return [str(e)]
	out = []
	if log['unreleased'] is not None and not log['unreleased']:
		out.append(f'"{UNRELEASED}" without entries: remove it until there is one')
	releases = log['releases']
	if not releases:
		return out + ['the changelog has no release']
	for entry in (log['unreleased'] or []) + [e for r in releases for e in r['entries']]:
		if not entry['text'].endswith(('.', '.)')):
			out.append(f'README.md:{entry["line"]}: an entry ends with a period')
	for release in releases:
		if not release['entries']:
			out.append(f'README.md:{release["line"]}: the release {release["tag"]} has no entries')
	for new, old in zip(releases, releases[1:]):
		if version(new['tag']) <= version(old['tag']):
			out.append(f'README.md:{new["line"]}: {new["tag"]} is above {old["tag"]}: newest first')
		if new['date'] < old['date']:
			out.append(f'README.md:{new["line"]}: {new["tag"]} has an older date than {old["tag"]}')
	known_tags = tags() if known_tags is None else known_tags
	if known_tags:
		listed = {version(r['tag']) for r in releases}
		missing = sorted((t for t in known_tags if version(t) not in listed), key=version)
		if missing:
			out.append(f'the releases {", ".join(missing)} (tags of git) have no heading in the changelog: a release '
			           f'renames "{UNRELEASED}" to its heading')
		tagged = {version(t) for t in known_tags}
		untagged = [r['tag'] for r in releases if version(r['tag']) not in tagged]
		if untagged:
			t = untagged[0]
			out.append(f'the releases {", ".join(untagged)} of the changelog have no tag of git, so no GitHub release '
			           f'nor version in Zenodo: push the tag (git push --follow-tags origin master) or, if it is not '
			           f'on the commit of its heading (e.g. after a rebase), make it there (git tag -d {t}; '
			           f'git tag -a {t} -m "eHabitat+ {t}" <commit>) and push it')
	return out


def entry_text(entry):
	"""An entry as a line of Markdown: '- **Desktop app** (fix): ...'."""
	return f'- **{", ".join(entry["areas"])}**{" (fix)" if entry["fix"] else ""}: {entry["text"]}'


def notes(text, tag):
	"""The notes of the GitHub release of TAG, in Markdown: its entries."""
	for release in parse(text)['releases']:
		if version(release['tag']) == version(tag):
			lines = [f'eHabitat+ {release["tag"]} ({release["date"]}).', '']
			lines += [entry_text(e) for e in release['entries']]
			lines += ['', f'All the versions: [changelog]({REPO}#changelog).']
			return '\n'.join(lines)
	raise ValueError(f'the changelog of README.md has no release {tag}')


def unrecorded(old, new):
	"""Why the commits OLD..NEW need an entry in the changelog and have none (the files of eHabitat+ that changed,
	with no line added to the section "## Changelog" of README.md), or None. Nothing to compare (a new branch, a
	commit that is not here): None."""
	if not old or set(old) == {'0'} or git('cat-file', '-e', f'{old}^{{commit}}') is None:
		return None
	files = (git('diff', '--name-only', old, new) or '').split()
	code = [f for f in files if CODE.match(f)]
	if not code:
		return None
	diff = git('diff', '-U0', old, new, '--', 'README.md') or ''
	added = [line[1:] for line in diff.splitlines() if line.startswith('+') and not line.startswith('+++')]
	if any(ENTRY.match(line.rstrip()) for line in added):
		return None
	shown = ', '.join(code[:5]) + (f' and {len(code) - 5} more' if len(code) > 5 else '')
	return (f'these commits change {shown} without an entry in the changelog of README.md: if users notice the '
	        f'change, add it under "{UNRELEASED}"')


def newest(text, known_tags=()):
	"""The newest version among the releases of the changelog and KNOWN_TAGS, as a tuple ((0, 0, 0) if none)."""
	found = [version(r['tag']) for r in parse(text)['releases']] + [version(t) for t in known_tags]
	return max(found, default=(0, 0, 0))


def next_version(text, known_tags=()):
	"""The release after the newest one: its last number plus one (1.6 -> 1.6.1, 1.6.1 -> 1.6.2)."""
	major, minor, patch = newest(text, known_tags)
	return f'{major}.{minor}.{patch + 1}'


def released(text, tag, date, known_tags=()):
	"""README.md with "### Unreleased" renamed to the heading of the release TAG of DATE (YYYY-MM-DD). Raises ValueError
	if there is nothing to release, TAG is not a version newer than the newest release or tag, or the changelog would
	have a problem."""
	if not re.fullmatch(r'(0|[1-9]\d*)\.(0|[1-9]\d*)(\.(0|[1-9]\d*))?', tag):
		raise ValueError(f'{tag} is not a version X.Y or X.Y.Z (numbers without leading zeros)')
	if not parse(text)['unreleased']:
		raise ValueError(f'nothing to release: "{UNRELEASED}" has no entries')
	last = newest(text, known_tags)
	if version(tag) <= last:
		raise ValueError(f'{tag} is not newer than the last release, {".".join(map(str, last))}')
	n = next(n for n, line in section(text) if line.strip() == UNRELEASED)
	lines = text.splitlines(keepends=True)
	end = lines[n - 1][len(lines[n - 1].rstrip('\r\n')):]
	lines[n - 1] = f'### [{tag}]({REPO}/releases/tag/{tag}) · {date}{end}'
	new = ''.join(lines)
	found = problems(new, [*known_tags, tag])
	if found:
		raise ValueError('; '.join(found))
	return new


def release(tag=None, trailers=()):
	"""The release of "### Unreleased" (see the top): renames it in README.md, commits README.md alone and tags that
	commit (annotated, so that git push --follow-tags pushes it). Returns the tag. Raises ValueError, and leaves
	README.md as it was, if it cannot be made."""
	if (git('rev-parse', '--abbrev-ref', 'HEAD') or '').strip() != 'master':
		raise ValueError('a release is made on master')
	if git('diff', '--quiet', 'HEAD', '--', 'README.md') is None:
		raise ValueError('README.md has changes that are not committed: commit them (with their entries) first')
	text, known = readme(), tags() or []
	tag = tag or next_version(text, known)
	if git('rev-parse', '-q', '--verify', f'refs/tags/{tag}') is not None:
		raise ValueError(f'the tag {tag} exists already')
	new = released(text, tag, datetime.datetime.now(datetime.timezone.utc).date().isoformat(), known)
	write_readme(new)
	extra = [a for t in trailers for a in ('--trailer', t)]
	err = git_error('commit', '-q', '-m', f'Release {tag}', *extra, '--', 'README.md')
	if err:
		write_readme(text)
		raise ValueError(f'git commit of the release {tag} failed ({err}): README.md is as it was')
	err = git_error('tag', '-a', tag, '-m', f'eHabitat+ {tag}')
	if err:
		raise ValueError(f'the commit "Release {tag}" is made, but git tag failed ({err}): make the tag with '
		                 f'git tag -a {tag} -m "eHabitat+ {tag}" (annotated, so that git push --follow-tags pushes it)')
	return tag


def git_error(*args):
	"""Runs git in the repository, without a time limit (a commit may wait for a passphrase): its error, or None if it
	worked."""
	try:
		r = subprocess.run(['git', '-C', ROOT, *args], capture_output=True, text=True)
	except OSError as e:
		return str(e)
	return None if r.returncode == 0 else (r.stderr.strip() or r.stdout.strip() or f'exit status {r.returncode}')


def pending(text):
	"""Why master has changes for users without their release (entries under "### Unreleased"), or None."""
	n = len(parse(text)['unreleased'] or [])
	if not n:
		return None
	return (f'"{UNRELEASED}" has {n} {"entry" if n == 1 else "entries"} pushed without a release: make it with python '
	        f'.github/changelog.py release and push it with git push --atomic --follow-tags origin master')


def readme():
	with open(os.path.join(ROOT, 'README.md'), encoding='utf-8', newline='') as f:
		return f.read()


def write_readme(text):
	with open(os.path.join(ROOT, 'README.md'), 'w', encoding='utf-8', newline='') as f:
		f.write(text)


def main(argv):
	usage = __doc__.split('\n\n')[1]
	if argv and argv[0] == 'release':
		args, trailers = list(argv[1:]), []
		while '--trailer' in args:
			i = args.index('--trailer')
			if i + 1 >= len(args):
				sys.exit(usage)
			trailers.append(args[i + 1])
			del args[i:i + 2]
		if len(args) > 1:
			sys.exit(usage)
		try:
			tag = release(args[0] if args else None, trailers)
		except ValueError as e:
			sys.exit(str(e))
		print(f'Release {tag} committed and tagged: push it with git push --atomic --follow-tags origin master')
		return
	cmds = {'check': 1, 'notes': 2, 'unrecorded': 3, 'pending': 1}
	if not argv or argv[0] not in cmds or len(argv) != cmds[argv[0]]:
		sys.exit(usage)
	try:
		if argv[0] == 'check':
			found = problems(readme())
			for p in found:
				print(p)
			sys.exit(1 if found else 0)
		if argv[0] == 'notes':
			print(notes(readme(), argv[1]))
		elif argv[0] == 'pending':
			why = pending(readme())
			if why:
				print(why)
		else:
			why = unrecorded(argv[1], argv[2])
			if why:
				print(why)
	except ValueError as e:
		sys.exit(str(e))


if __name__ == '__main__':
	main(sys.argv[1:])
