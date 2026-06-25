# Releasing a new version

This is the maintainer checklist for cutting a new pyMARS release. Releases are
published to **PyPI** (as `nrg-pymars`) and **conda-forge**; the import package
and command-line tool are named `pymars`.

## Prerequisites (one-time setup)

- **PyPI:** OIDC trusted publishing is configured for the `nrg-pymars` project,
  tied to this repository's `publish.yml` workflow and the `pypi` environment.
- **conda-forge:** the `nrg-pymars-feedstock` exists (created automatically after
  the `conda-forge/staged-recipes` PR is merged).

## Release steps

1. Make sure `main` is green (CI passing).

2. Bump the version in `src/pymars/_version.py`. This is the single source of
   truth — `pyproject.toml`, the conda recipe, and the docs all read from it.

3. Update `CHANGELOG.md`: move the entries under `## [Unreleased]` into a new
   `## [X.Y.Z] - YYYY-MM-DD` section, and add the matching compare link at the
   bottom of the file.

4. Update `CITATION.cff`: set `version` and `date-released`.

5. Open a pull request with the above changes, wait for CI to pass, and merge it
   (branch protection requires a PR to land on `main`).

6. Tag the merge commit on `main` and push the tag:

   ```bash
   git switch main
   git pull upstream main
   git tag -a vX.Y.Z -m "vX.Y.Z"
   git push upstream vX.Y.Z
   ```

   Pushing the tag triggers `.github/workflows/publish.yml`, which runs the CI
   test suite first and then, only if it passes, builds and publishes to PyPI.

7. **conda-forge** (mostly automatic): within a few hours of the PyPI release,
   the conda-forge autotick bot opens a version-bump pull request on the
   `nrg-pymars-feedstock` repository, updating the version and sha256. Review and
   merge it; conda-forge CI then builds and publishes to the `conda-forge`
   channel. You only need to edit the recipe yourself if the dependencies or the
   minimum Python version changed.

8. (Optional) Create a GitHub Release for the tag, using the new `CHANGELOG.md`
   section as the release notes.

9. Verify the release:

   ```bash
   pip install nrg-pymars==X.Y.Z
   # once the conda-forge build completes:
   conda install -c conda-forge nrg-pymars
   ```

## TODO: clean up once the conda-forge package is live

The self-hosted Anaconda.org build is redundant once `nrg-pymars` is on
conda-forge. After the feedstock is publishing, do the following:

- [ ] Remove the `conda` job from `.github/workflows/publish.yml` (the
      Anaconda.org build + upload to the `niemeyer-research-group` channel).
- [ ] Remove `conda.recipe/` (the self-hosted recipe) if nothing else uses it.
- [ ] Update the conda install instructions in `README.md` and
      `docs/installation.rst` to `conda install -c conda-forge nrg-pymars`.
- [ ] Remove the `ANACONDA_API_TOKEN` repository secret once nothing references it.
- [ ] Update this file to drop the Anaconda.org references above.
