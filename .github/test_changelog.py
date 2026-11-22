"""Tests of .github/changelog.py: the next version, the renaming of "### Unreleased" (released), the warning of entries
pushed without a release (pending) and the release itself (commit and annotated tag) in a temporary git repository.
Only the standard library, like changelog.py; it runs in .github/workflows/changelog.yml.

    python .github/test_changelog.py
"""

import datetime
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import changelog as cl  # noqa: E402

REPO = cl.REPO
README = f"""# eHabitat+

## Changelog

What each release brings.

### Unreleased

- **HFT model**: Something new.
- **Desktop app** (fix): Something fixed.

### [1.5]({REPO}/releases/tag/1.5) · 2026-11-14

- **Desktop app**: Older thing.

### [1.4]({REPO}/releases/tag/1.4) · 2026-10-27

- **HFT model**: Oldest thing.
"""
NO_UNRELEASED = README.replace('### Unreleased\n\n- **HFT model**: Something new.\n- **Desktop app** (fix): Something '
                               'fixed.\n\n', '')
TODAY = '2026-11-22'
TAGS = ['1.4', '1.5']


class Versions(unittest.TestCase):
	def test_next(self):
		self.assertEqual(cl.next_version(README), '1.5.1')
		self.assertEqual(cl.next_version(README, ['1.4', '1.5', '1.6']), '1.6.1')
		self.assertEqual(cl.next_version(README, ['1.6.1']), '1.6.2')
		self.assertEqual(cl.next_version(README, ['v1.6.9']), '1.6.10')


class Released(unittest.TestCase):
	def test_renames_unreleased(self):
		new = cl.released(README, '1.6', TODAY, ['1.4', '1.5'])
		self.assertNotIn('### Unreleased', new)
		self.assertIn(f'### [1.6]({REPO}/releases/tag/1.6) · {TODAY}\n\n- **HFT model**: Something new.', new)
		self.assertEqual(cl.problems(new, ['1.4', '1.5', '1.6']), [])
		self.assertEqual(new.replace(f'### [1.6]({REPO}/releases/tag/1.6) · {TODAY}', '### Unreleased'), README)
		notes = cl.notes(new, '1.6')
		self.assertTrue(notes.startswith(f'eHabitat+ 1.6 ({TODAY}).'))
		self.assertIn('- **Desktop app** (fix): Something fixed.', notes)

	def test_keeps_crlf(self):
		crlf = README.replace('\n', '\r\n')
		new = cl.released(crlf, '1.6', TODAY, TAGS)
		self.assertIn(f'### [1.6]({REPO}/releases/tag/1.6) · {TODAY}\r\n', new)
		self.assertNotIn('\n\n', new.replace('\r\n', '#'))

	def test_refuses(self):
		for tag, known, text, why in [
				('1.6', TAGS, NO_UNRELEASED, 'nothing to release'),
				('1.5', TAGS, README, 'not newer'),                       # not newer than the last heading
				('1.4.9', TAGS, README, 'not newer'),
				('1.6', [*TAGS, '1.6'], README, 'not newer'),             # not newer than a tag
				('1.6', TAGS, README.replace('2026-11-14', '2026-12-01'), 'older date'),
				('1.6', ['1.4'], README, 'have no tag of git'),           # a heading without its tag (1.5)
				('v1.6', TAGS, README, 'not a version'), ('1.6.0.1', TAGS, README, 'not a version'),
				('1.6-beta', TAGS, README, 'not a version'), ('1.07', TAGS, README, 'leading zeros'),
				('01.7', TAGS, README, 'leading zeros')]:
			with self.subTest(tag=tag, known=known):
				with self.assertRaises(ValueError) as e:
					cl.released(text, tag, TODAY, known)
				self.assertIn(why, str(e.exception))

	def test_empty_unreleased(self):
		with self.assertRaises(ValueError):
			cl.released(README.replace('- **HFT model**: Something new.\n- **Desktop app** (fix): Something fixed.\n',
			                           ''), '1.6', TODAY, TAGS)


class Pending(unittest.TestCase):
	def test_pending(self):
		self.assertIn('2 entries pushed without a release', cl.pending(README))
		self.assertIsNone(cl.pending(NO_UNRELEASED))
		self.assertIsNone(cl.pending(cl.released(README, '1.6', TODAY, TAGS)))


class Release(unittest.TestCase):
	"""changelog.py release in a temporary repository with its own copy of changelog.py."""

	def setUp(self):
		self.dir = tempfile.mkdtemp()
		os.makedirs(os.path.join(self.dir, '.github'))
		shutil.copy(os.path.join(HERE, 'changelog.py'), os.path.join(self.dir, '.github'))
		self.write(README)
		self.git('init', '-q', '-b', 'master')
		self.git('config', 'user.name', 'Test')
		self.git('config', 'user.email', 'test@example.org')
		self.git('config', 'commit.gpgsign', 'false')
		self.git('config', 'tag.gpgsign', 'false')
		self.git('add', '.')
		self.git('commit', '-q', '-m', 'change')
		self.git('tag', '1.4')
		self.git('tag', '1.5')

	def tearDown(self):
		shutil.rmtree(self.dir)

	def git(self, *args):
		return subprocess.run(['git', '-C', self.dir, *args], capture_output=True, text=True, check=True).stdout

	def write(self, text):
		with open(os.path.join(self.dir, 'README.md'), 'w', encoding='utf-8', newline='') as f:
			f.write(text)

	def read(self):
		with open(os.path.join(self.dir, 'README.md'), encoding='utf-8', newline='') as f:
			return f.read()

	def run_cl(self, *args):
		return subprocess.run([sys.executable, os.path.join(self.dir, '.github', 'changelog.py'), *args],
		                      capture_output=True, text=True)

	def test_release(self):
		r = self.run_cl('release', '1.6', '--trailer', 'Co-Authored-By: Someone <someone@example.org>')
		self.assertEqual(r.returncode, 0, r.stderr)
		today = datetime.datetime.now(datetime.timezone.utc).date().isoformat()
		self.assertIn(f'### [1.6]({REPO}/releases/tag/1.6) · {today}', self.read())
		self.assertEqual(self.git('log', '-1', '--format=%B').strip(),
		                 'Release 1.6\n\nCo-Authored-By: Someone <someone@example.org>')
		self.assertEqual(self.git('diff', '--name-only', 'HEAD~1', 'HEAD').split(), ['README.md'])
		self.assertEqual(self.git('cat-file', '-t', '1.6').strip(), 'tag')               # annotated
		self.assertEqual(self.git('rev-parse', '1.6^{commit}'), self.git('rev-parse', 'HEAD'))
		self.assertEqual(self.git('status', '--porcelain'), '')
		self.assertEqual(self.run_cl('check').returncode, 0)
		# nothing more to release: it fails and leaves README.md as it is
		before = self.read()
		r = self.run_cl('release')
		self.assertNotEqual(r.returncode, 0)
		self.assertIn('nothing to release', r.stderr)
		self.assertEqual(self.read(), before)
		# a new entry: the next release is 1.6.1, then 1.6.2
		for nxt in ('1.6.1', '1.6.2'):
			self.write(self.read().replace('## Changelog\n\nWhat each release brings.\n\n',
			                               f'## Changelog\n\nWhat each release brings.\n\n### Unreleased\n\n'
			                               f'- **Similarity**: New for {nxt}.\n\n'))
			self.git('commit', '-q', '-am', f'change for {nxt}')
			r = self.run_cl('release')
			self.assertEqual(r.returncode, 0, r.stderr)
			self.assertIn(f'Release {nxt} committed and tagged', r.stdout)
			self.assertEqual(self.git('describe', '--tags', '--exact-match').strip(), nxt)

	def test_refuses_uncommitted_readme(self):
		head = self.git('rev-parse', 'HEAD')
		edited = README.replace('Something new.', 'Something newer.')
		self.write(edited)
		for staged in (False, True):
			if staged:
				self.git('add', 'README.md')
			r = self.run_cl('release')
			self.assertIn('not committed', r.stderr)
			self.assertEqual(self.read(), edited)
			self.assertEqual(self.git('rev-parse', 'HEAD'), head)
			self.assertEqual(self.git('tag').split(), ['1.4', '1.5'])

	def test_failed_commit_leaves_readme(self):
		hook = os.path.join(self.dir, '.git', 'hooks', 'pre-commit')
		with open(hook, 'w') as f:
			f.write('#!/bin/sh\necho the hook says no >&2\nexit 1\n')
		os.chmod(hook, 0o755)
		head = self.git('rev-parse', 'HEAD')
		r = self.run_cl('release', '1.6')
		self.assertIn('the hook says no', r.stderr)
		self.assertIn('README.md is as it was', r.stderr)
		self.assertEqual(self.read(), README)
		self.assertEqual(self.git('rev-parse', 'HEAD'), head)
		self.assertEqual(self.git('status', '--porcelain'), '')
		self.assertEqual(self.git('tag').split(), ['1.4', '1.5'])

	def test_failed_tag_says_how(self):
		if not os.path.isdir(os.path.join(self.dir, '.git', 'refs', 'tags')):
			self.skipTest('git stores its refs in a reftable: no lock file to make git tag fail')
		with open(os.path.join(self.dir, '.git', 'refs', 'tags', '1.6.lock'), 'w'):
			pass  # git tag 1.6 cannot lock its ref
		r = self.run_cl('release', '1.6')
		self.assertNotEqual(r.returncode, 0)
		self.assertIn('git tag -a 1.6 -m "eHabitat+ 1.6"', r.stderr)
		self.assertEqual(self.git('log', '-1', '--format=%s').strip(), 'Release 1.6')
		# the heading without its tag is a problem of check
		r = self.run_cl('check')
		self.assertNotEqual(r.returncode, 0)
		self.assertIn('1.6 of the changelog have no tag of git', r.stdout)

	def test_refuses_other_branch(self):
		self.git('checkout', '-q', '-b', 'other')
		r = self.run_cl('release')
		self.assertIn('on master', r.stderr)
		self.assertEqual(self.read(), README)

	def test_refuses_existing_tag(self):
		self.git('tag', '1.5.1')
		r = self.run_cl('release', '1.5.1')
		self.assertNotEqual(r.returncode, 0)
		self.assertEqual(self.read(), README)
		self.assertEqual(self.git('log', '-1', '--format=%s').strip(), 'change')

	def test_usage(self):
		for args in (['release', '1.6', '1.7'], ['release', '--trailer'], ['pending', 'x']):
			with self.subTest(args=args):
				r = self.run_cl(*args)
				self.assertNotEqual(r.returncode, 0)
				self.assertIn('changelog.py check', r.stderr)


if __name__ == '__main__':
	unittest.main(verbosity=1)
