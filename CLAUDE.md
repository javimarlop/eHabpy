# eHabitat+ (eHabpy)

## Every change on its own branch, merged into master before the push

The user may make changes in parallel (several sessions, or the user and Claude at once), so nothing reaches GitHub
before it is merged into `master` here (the user's rule, 2026-11-22):

- Every change is made on its own branch from the latest `origin/master`, in its own folder (a worktree: a folder has
  one branch checked out, so two branches in the same folder would switch each other's files):
  `git fetch origin && git worktree add ../eHabpy-<branch> -b <branch> origin/master`, or a Claude Code session in its
  own worktree. Never commit a change on `master` itself.
- This folder (`~/Downloads/Claude/eHabpy`) stays on `master`: here branches are only merged, released and pushed.
  Never switch it to a branch, nor leave uncommitted changes in it.
- A branch is never pushed to GitHub, and there are no pull requests (the user works alone): its commits reach GitHub
  only through `master`.
- When the change is done, tested and, if it reaches users, reviewed, and the user says to push it: in this folder,
  `git pull --ff-only origin master`, then `git merge --no-ff <branch>`. Resolve the conflicts (two changes with
  entries at the top of `### Unreleased`: keep both), run the tests of what the merge touched and
  `python .github/changelog.py check`, then make its release and push (below). Then
  `git worktree remove ../eHabpy-<branch>` and `git branch -d <branch>`.
- To bring the changes merged by others into a branch before its merge: `git merge origin/master` in its worktree.

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
- **Every push to `master` with entries under `### Unreleased` is a release** (the user's decision, 2026-11-22): 1.6,
  then 1.6.1, 1.6.2... (the last number plus one; a new minor version, e.g. `release 1.7`, only when the user asks). A
  push with nothing for users (no entries) makes none. After the merge of its branch into `master`, in this folder:
  `python .github/changelog.py release` renames `### Unreleased` to
  `### [X.Y.Z](https://github.com/javimarlop/eHabpy/releases/tag/X.Y.Z) · YYYY-MM-DD` (today in UTC), commits
  `README.md` alone (`Release X.Y.Z`; `--trailer "..."` adds the attribution line) and makes an annotated tag; then
  `git push --atomic --follow-tags origin master`, one release per push. The tag makes `.github/workflows/release.yml`
  create the GitHub release `eHabitat+ X.Y.Z`, whose notes are the output of `python .github/changelog.py notes X.Y.Z`,
  and Zenodo archives it with its own DOI. A version of Zenodo cannot be deleted: check the entries before the push.
  If the push is rejected (master moved), undo the release (`git tag -d X.Y.Z; git reset --keep HEAD~1`), pull and
  make it again; `check` fails on a heading of a release without its tag. By hand, Actions > Release > Run workflow
  with a tag only checks and shows its notes, or creates a missing release with "publish".
- Zenodo archives every GitHub release (concept DOI 10.5281/zenodo.595844). `.zenodo.json` keeps the same title,
  author (with ORCID) and licence (CC BY-SA 3.0, that of `LICENSE`) in every release; the description of each one is
  its release notes. Do not add `description` or `related_identifiers` to `.zenodo.json`: they would replace the notes
  and the link to the code of each release.
- `python .github/changelog.py check` checks the format; the workflow `.github/workflows/changelog.yml` runs it and
  `.github/test_changelog.py` on every push, warns when the commits change the app, the scripts or the environment
  without an entry, and when a push to `master` leaves entries under `### Unreleased` without their release.
