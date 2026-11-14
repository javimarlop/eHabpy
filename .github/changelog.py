"""The changelog of eHabitat+: the section "## Changelog" at the end of README.md, what each release brings, newest
first (as that of MonitorEO, https://github.com/monitoreo-obsnev/monitoreo#changelog).

    python .github/changelog.py check                 # its format, and that it has the releases of the tags
    python .github/changelog.py notes 1.5             # the entries of a release in Markdown (its GitHub release notes)
    python .github/changelog.py unrecorded OLD NEW    # whether the commits OLD..NEW change eHabitat+ without an entry

Each release is a heading with its tag and the date of its GitHub release, followed by its entries, one per line, each
with the parts of eHabitat+ it reaches (AREAS) and "(fix)" for the fixes:

    ### [1.4](https://github.com/javimarlop/eHabpy/releases/tag/1.4) · 2026-10-27

    - **HFT model**: New CLASS_METHOD option in config.py ...
    - **Similarity** (fix): The similarity indices are right when a mask is used.

Above the releases, "### Unreleased" holds the entries of the changes pushed to master since the last release (the
desktop app has them already: it updates itself from master every time it opens); a release renames it. Every change
that users notice adds its entry there, in the same commit. Only the standard library: it runs in GitHub Actions
(.github/workflows/changelog.yml).
"""

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


def readme():
	with open(os.path.join(ROOT, 'README.md'), encoding='utf-8') as f:
		return f.read()


def main(argv):
	cmds = {'check': 1, 'notes': 2, 'unrecorded': 3}
	if not argv or argv[0] not in cmds or len(argv) != cmds[argv[0]]:
		sys.exit(__doc__.split('\n\n')[1])
	try:
		if argv[0] == 'check':
			found = problems(readme())
			for p in found:
				print(p)
			sys.exit(1 if found else 0)
		if argv[0] == 'notes':
			print(notes(readme(), argv[1]))
		else:
			why = unrecorded(argv[1], argv[2])
			if why:
				print(why)
	except ValueError as e:
		sys.exit(str(e))


if __name__ == '__main__':
	main(sys.argv[1:])
