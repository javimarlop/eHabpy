"""Tests of .github/wiki_links.py: the anchors GitHub gives to the headings and the broken links found in a small wiki
written in a temporary folder. Only the standard library, like wiki_links.py; it runs in .github/workflows/wiki.yml.

    python .github/test_wiki_links.py
"""

import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wiki_links as wl  # noqa: E402

PAGES = {
	'Home.md': """# eHabitat+
<p><img src="images/logo.png"></p>

## Contents
- [Tutorial](HFT-tutorial#step-6-run-the-model-run_hftpy), [plots](HFT-tutorial#variance-explained-r2_id_hftspng)
- [[HFT tutorial]], [[The steps|HFT-tutorial]], [this page](#contents), [the home](https://github.com/javimarlop/eHabpy/wiki)
- [Paper](https://doi.org/10.1/x), [code](../blob/master/README.md)
- Broken: [a](Missing-page), [b](HFT-tutorial#no-such-step), [c](#nowhere), [d](images/none.png), [e](hft-tutorial)
- `[in code](Missing-in-code)`
""",
	'HFT-tutorial.md': """## Step 6. Run the model (run_hft.py)
## Variance explained (`r2_<id>_hfts.png`)
## Troubleshooting
## Troubleshooting
<a name="by-hand"></a>
```bash
# Not a heading
[not a link](Missing-in-a-block)
```
[the second](#troubleshooting-1), [by hand](#by-hand)
""",
	'_Sidebar.md': """<details>
<summary><a href="https://github.com/javimarlop/eHabpy/wiki/HFT-tutorial">HFT tutorial</a></summary>

- [Steps](HFT-tutorial#troubleshooting)
- <a href="https://github.com/javimarlop/eHabpy/wiki/HFT-tutorial#gone">gone</a>

</details>
""",
}


class Slug(unittest.TestCase):

	def test_headings(self):
		self.assertEqual(wl.slug('Step 6. Run the model (run_hft.py)'), 'step-6-run-the-model-run_hftpy')
		self.assertEqual(wl.slug('Variance explained by the HFTs (`r2_<id>_hfts.png`)'),
		                 'variance-explained-by-the-hfts-r2_id_hftspng')
		self.assertEqual(wl.slug('A12. Similarity of the HFTs, in [the ecoregion](x)'),
		                 'a12-similarity-of-the-hfts-in-the-ecoregion')
		self.assertEqual(wl.slug('**Bold** <b>tag</b> χ² ñ'), 'bold-tag-χ²-ñ')

	def test_repeated_and_code(self):
		found = wl.anchors(PAGES['HFT-tutorial.md'])
		self.assertIn('troubleshooting', found)
		self.assertIn('troubleshooting-1', found)
		self.assertIn('by-hand', found)
		self.assertNotIn('not-a-heading', found)


class Links(unittest.TestCase):

	def setUp(self):
		self.dir = tempfile.mkdtemp()
		os.makedirs(os.path.join(self.dir, 'images'))
		open(os.path.join(self.dir, 'images', 'logo.png'), 'w').close()
		for name, text in PAGES.items():
			with open(os.path.join(self.dir, name), 'w', encoding='utf-8') as f:
				f.write(text)
		self.wiki = wl.Wiki(self.dir)

	def tearDown(self):
		shutil.rmtree(self.dir)

	def broken(self, page):
		return [(where.split(':')[1], link) for where, link, _ in
		        wl.check_page(self.wiki, page[:-3], PAGES[page], page)]

	def test_pages(self):
		self.assertEqual(self.broken('Home.md'), [('8', 'Missing-page'), ('8', 'HFT-tutorial#no-such-step'),
		                                          ('8', '#nowhere'), ('8', 'images/none.png'), ('8', 'hft-tutorial')])
		self.assertEqual(self.broken('HFT-tutorial.md'), [])
		self.assertEqual(self.broken('_Sidebar.md'),
		                 [('5', 'https://github.com/javimarlop/eHabpy/wiki/HFT-tutorial#gone')])

	def test_why(self):
		self.assertEqual(self.wiki.problem('hft-tutorial'), 'no page hft-tutorial (there is HFT-tutorial)')
		self.assertEqual(self.wiki.problem('#nowhere', 'Home'), 'no section #nowhere in Home')

	def test_files_of_the_repository(self):
		path = os.path.join(self.dir, 'README.txt')
		with open(path, 'w') as f:
			f.write('See https://github.com/javimarlop/eHabpy/wiki/HFT-tutorial#troubleshooting.\n'
			        '(https://github.com/javimarlop/eHabpy/wiki#contents) and '
			        'https://github.com/javimarlop/eHabpy/wiki/Old-page\n')
		self.assertEqual([(w.split(':')[-1], u) for w, u, _ in wl.check_file(self.wiki, path)],
		                 [('2', 'https://github.com/javimarlop/eHabpy/wiki/Old-page')])

	def test_main(self):
		self.assertEqual(wl.main([self.dir]), 1)
		os.remove(os.path.join(self.dir, 'Home.md'))
		os.remove(os.path.join(self.dir, '_Sidebar.md'))
		self.assertEqual(wl.main([self.dir]), 0)


if __name__ == '__main__':
	unittest.main()
