# Contributing

We welcome contributions in the form of bug reports, bug fixes, improvements to the documentation, ideas for enhancements, or the enhancements themselves!

You can find a [list of current issues](https://github.com/Niemeyer-Research-Group/pyMARS/issues) in the project's GitHub repository. Feel free to tackle any existing bugs or enhancement ideas by submitting a [pull request](https://github.com/Niemeyer-Research-Group/pyMARS/pulls). Some issues are marked as `beginner-friendly`. These issues are a great place to start working with pyMARS, if you're new here.

## Bug Reports

 * Please include a short (but detailed) Python snippet or explanation for reproducing the problem. Attach or include a link to any input files that will be needed to reproduce the error.
 * Explain the behavior you expected, and how what you got differed.
 * Include the full text of any error messages that are printed on the screen.

## Development setup

We recommend working inside an isolated environment (see the [installation guide](https://Niemeyer-Research-Group.github.io/pyMARS/installation.html)). Clone the repository, install pyMARS with its development dependencies, and set up the pre-commit hooks:

    git clone https://github.com/Niemeyer-Research-Group/pyMARS.git
    cd pyMARS
    pip install -e .[dev]
    pre-commit install

## Pull Requests

 * If you're unfamiliar with Pull Requests, please take a look at the [GitHub documentation for them](https://help.github.com/articles/proposing-changes-to-a-project-with-pull-requests/).
 * Start by creating a new branch from the latest commit on [`main`](https://github.com/Niemeyer-Research-Group/pyMARS/tree/main).
 * **Make sure the test suite passes** and that test coverage doesn't go down. From the top-level directory, run `pytest --cov=pymars`.
 * *Always* add tests and docs for your code.
 * Code style is enforced with [Black](https://black.readthedocs.io/) and [Ruff](https://docs.astral.sh/ruff/) via [pre-commit](https://pre-commit.com/). Run `pre-commit run --all-files` (or install the hooks as shown above) before pushing.
 * Docstrings are required and should follow the [NumPy style](https://numpydoc.readthedocs.io/en/latest/format.html).
 * Add an entry describing your changes to the [`CHANGELOG`](https://github.com/Niemeyer-Research-Group/pyMARS/blob/main/CHANGELOG.md), under the `Unreleased` section.
 * The use of emoji in Pull Request titles is encouraged, with the format ":emoji: Commit summary". See [this list of suggested emoji](https://github.com/slashsBin/styleguide-git-commit-message#suggested-emojis).
 * Please reference relevant GitHub issues in your commit messages using `#123`.
 * The copyright policy is detailed in the [`LICENSE`](https://github.com/Niemeyer-Research-Group/pyMARS/blob/main/LICENSE).

## Releasing

If you're a maintainer cutting a new release, see [`RELEASING.md`](RELEASING.md) for the step-by-step checklist (version bump, changelog, tagging, and publishing to PyPI and conda-forge).

## Meta

Thanks to the useful [contributing guide of pyrk](https://github.com/pyrk/pyrk/blob/master/CONTRIBUTING.md), which served as an inspiration and starting point for this guide.
