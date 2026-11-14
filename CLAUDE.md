# eHabitat+ (eHabpy)

## Users see every change in the changelog

The changelog is the section `## Changelog` at the end of `README.md`, as that of MonitorEO
(https://github.com/monitoreo-obsnev/monitoreo#changelog). There is no other (the old `changelog.txt` is in its
history).

- A change that users notice adds its entry, in English, at the top of `### Unreleased` (create the heading above the
  releases if it is not there), **in the same commit as the change**: one line each,
  `- **Desktop app, HFT model**: <what users get>.`, or `- **Similarity** (fix): <the error that no longer happens>.`
  Say what users get, not how it was done. The areas are those of `AREAS` in `.github/changelog.py`: eHabitat+,
  Desktop app, Installers, HFT model, Similarity, Earth Engine variables, Documentation.
- The desktop app updates itself from `master` every time it opens, so what is under `### Unreleased` already reaches
  users.
- A release (a tag and a GitHub release, archived in Zenodo) renames `### Unreleased` to
  `### [X.Y](https://github.com/javimarlop/eHabpy/releases/tag/X.Y) · YYYY-MM-DD`, and its notes are the output of
  `python .github/changelog.py notes X.Y`.
- Zenodo archives every GitHub release (concept DOI 10.5281/zenodo.595844). `.zenodo.json` keeps the same title,
  author (with ORCID) and licence (CC BY-SA 3.0, that of `LICENSE`) in every release; the description of each one is
  its release notes. Do not add `description` or `related_identifiers` to `.zenodo.json`: they would replace the notes
  and the link to the code of each release.
- `python .github/changelog.py check` checks the format; the workflow `.github/workflows/changelog.yml` runs it on every
  push and warns when the commits change the app, the scripts or the environment without an entry.
