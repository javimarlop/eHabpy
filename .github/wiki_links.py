"""The links of the wiki of eHabitat+ (https://github.com/monitoreo-obsnev/eHabpy/wiki): every link to a page, a section
(anchor) or an image of the wiki must lead somewhere, in the pages of the wiki (_Sidebar.md and _Footer.md included)
and in the files of the repository that link to the wiki (README.md, the comments of the scripts...).

    python .github/wiki_links.py WIKI [FILE ...]

WIKI is a copy of the wiki (git clone https://github.com/monitoreo-obsnev/eHabpy.wiki.git), FILE the files of the repository
whose links to https://github.com/monitoreo-obsnev/eHabpy/wiki are also checked. It prints one line per broken link
("file:line: link: why") and exits with 1 if there is any, 0 otherwise. Only the standard library; it runs in
.github/workflows/wiki.yml every time the wiki is edited.

A section is reached as GitHub names it: the text of its heading without its Markdown (links, code, emphasis) nor its
HTML, in lower case, without the characters that are not letters, digits, spaces, "-" or "_", and with "-" for every
space; a heading repeated in a page gets "-1", "-2"... Headings and links inside code (fenced blocks and `inline code`)
do not count. Anchors written by hand (<a name="x">, id="x") count too. Headings underlined with === or --- (setext)
are not read: the wiki uses only # headings.
"""

import os
import re
import sys
import urllib.parse

WIKI_URL = 'https://github.com/monitoreo-obsnev/eHabpy/wiki'
HEADING = re.compile(r'^ {0,3}(#{1,6})\s+(.*?)(?:\s+#+)?\s*$')
FENCE = re.compile(r'^ {0,3}(`{3,}|~{3,})')
INLINE_CODE = re.compile(r'(`+)(.+?)\1')
MD_LINK = re.compile(r'!?\[((?:[^\[\]]|\[[^\]]*\])*)\]\(\s*<?([^)\s>]*)>?(?:\s+["\'(][^)]*)?\)')
WIKI_LINK = re.compile(r'\[\[([^\]]+)\]\]')
HTML_LINK = re.compile(r'\b(?:href|src)\s*=\s*["\']([^"\']*)["\']', re.I)
HTML_ANCHOR = re.compile(r'<a\s[^>]*\bname\s*=\s*["\']([^"\']+)["\']|\bid\s*=\s*["\']([^"\']+)["\']', re.I)
URL_IN_TEXT = re.compile(re.escape(WIKI_URL) + r'[^\s)\]"\'<>`]*')


def slug(heading):
	"""The anchor that GitHub gives to a heading: "Step 6. Run the model (run_hft.py)" -> step-6-run-the-model-run_hftpy."""
	parts = INLINE_CODE.split(heading)                             # text, ticks, code, text, ticks, code... text
	text = ''
	for i, part in enumerate(parts):
		if i % 3 == 2:                                                # code: as it is (<id> in code is text)
			text += part
		elif i % 3 == 0:                                              # text: without its Markdown and HTML
			part = re.sub(r'!\[[^\]]*\]\([^)]*\)', '', part)            # images
			part = re.sub(r'\[([^\]]*)\]\([^)]*\)', r'\1', part)        # links: their text
			text += re.sub(r'<[^>]+>', '', part)                      # HTML tags
	text = text.strip().lower()
	text = re.sub(r'[^\w\- ]', '', text)                            # \w: letters (also accented), digits and _
	return text.replace(' ', '-')


def lines_outside_code(text):
	"""(number, line, line with its `inline code` blanked) of the lines that are not in a fenced code block."""
	fence = None
	for n, line in enumerate(text.splitlines(), 1):
		m = FENCE.match(line)
		if fence:
			if m and m.group(1)[0] == fence[0] and len(m.group(1)) >= len(fence) and not line[m.end():].strip():
				fence = None
			continue
		if m:
			fence = m.group(1)
			continue
		yield n, line, INLINE_CODE.sub(lambda c: ' ' * len(c.group(0)), line)


def anchors(text):
	"""The anchors of a page: those of its headings (repeated ones with -1, -2...) and those written in HTML."""
	seen, out = {}, set()
	for _, line, bare in lines_outside_code(text):
		m = HEADING.match(line)                                       # the heading with its code: it is in the anchor
		if m:
			base = slug(m.group(2))
			k = seen.get(base, 0)
			seen[base] = k + 1
			out.add(base if k == 0 else f'{base}-{k}')
		for a, b in HTML_ANCHOR.findall(bare):
			out.add(a or b)
	return out


class Wiki:
	"""The pages (name -> their anchors) and the files of a copy of the wiki."""

	def __init__(self, folder):
		self.folder = folder
		self.pages = {}
		for name in sorted(os.listdir(folder)):
			if name.endswith('.md'):
				with open(os.path.join(folder, name), encoding='utf-8') as f:
					self.pages[name[:-3]] = anchors(f.read())

	def problem(self, target, page=None):
		"""Why a link of the wiki (Page, Page#anchor, #anchor, images/x.png; page = where it is, for #anchor) leads
		nowhere, or None if it leads somewhere."""
		target = urllib.parse.unquote(target)
		name, _, anchor = target.partition('#')
		name = name.rstrip('/') or (page if page is not None else 'Home')
		if name not in self.pages:
			if '.' in os.path.basename(name):                         # a file of the wiki (images/logo.png)
				return None if os.path.isfile(os.path.join(self.folder, name)) else f'no file {name} in the wiki'
			near = [p for p in self.pages if p.lower() == name.lower()]
			return f'no page {name}' + (f' (there is {near[0]})' if near else '')
		if anchor and anchor.lower() not in self.pages[name]:
			return f'no section #{anchor} in {name}'
		return None


def wiki_target(url):
	"""The part of a link to the wiki after .../wiki/ ('' for its home), or None if it does not lead to the wiki."""
	if url == WIKI_URL or url.startswith(WIKI_URL + '#'):
		return 'Home' + url[len(WIKI_URL):]
	if url.startswith(WIKI_URL + '/'):
		return url[len(WIKI_URL) + 1:] or 'Home'
	return None


def links_of_page(text):
	"""(line, link) of the links of a wiki page: Markdown ([x](y), ![x](y)), [[Page]], [[Text|Page]] and HTML (href,
	src)."""
	for n, _, line in lines_outside_code(text):
		for _, target in MD_LINK.findall(line):
			yield n, target
		for inner in WIKI_LINK.findall(line):
			yield n, '[[' + inner + ']]'
		for target in HTML_LINK.findall(line):
			yield n, target


def check_page(wiki, page, text, where):
	"""The broken links of a wiki page: (where:line, link, why)."""
	out = []
	for n, link in links_of_page(text):
		if link.startswith('[['):
			target = link[2:-2].split('|')[-1].strip().replace(' ', '-')
		elif re.match(r'^[a-z][a-z0-9+.-]*:', link, re.I) or link.startswith('//'):
			target = wiki_target(link)
			if target is None:                                        # another site: not checked
				continue
		elif link.startswith('/') or link.startswith('../') or not link:
			continue                                                  # other places of github.com, empty
		else:
			target = link
		why = wiki.problem(target, page)
		if why:
			out.append((f'{where}:{n}', link, why))
	return out


def check_file(wiki, path):
	"""The broken links to the wiki (https://github.com/monitoreo-obsnev/eHabpy/wiki/...) of a file of the repository."""
	out = []
	with open(path, encoding='utf-8', errors='replace') as f:
		for n, line in enumerate(f, 1):
			for url in URL_IN_TEXT.findall(line):
				url = url.rstrip('.,;:')
				why = wiki.problem(wiki_target(url))
				if why:
					out.append((f'{path}:{n}', url, why))
	return out


def main(argv):
	if not argv or argv[0] in ('-h', '--help'):
		print(__doc__)
		return 0 if argv else 2
	wiki = Wiki(argv[0])
	if not wiki.pages:
		print(f'ERROR: no pages (.md) in {argv[0]}')
		return 2
	problems = []
	for page in wiki.pages:
		with open(os.path.join(argv[0], page + '.md'), encoding='utf-8') as f:
			problems += check_page(wiki, page, f.read(), f'wiki/{page}.md')
	for path in argv[1:]:
		problems += check_file(wiki, path)
	for where, link, why in problems:
		print(f'{where}: {link}: {why}')
	print(f'{len(problems)} broken link(s) in {len(wiki.pages)} wiki pages and {len(argv) - 1} other files')
	return 1 if problems else 0


if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))
