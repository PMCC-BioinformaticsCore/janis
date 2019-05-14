# Releasing Janis

[![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis) [![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines)

Releasing janis is straight forward. Decide on a logical set of changes to include, this will 

### Versioning

Janis follows the  [`SemVer`](https://semver.org) versioning system:

> 1.  MAJOR version when you make incompatible API changes,
> 2.  MINOR version when you add functionality in a backwards-compatible manner, and
> 3.  PATCH version when you make backwards-compatible bug fixes.

Before a new release, you should update the version tag in `janis/__init__.py` so the produced python package is 

### Tagging and Building

> Before you tag, make sure you've incremented the `__version__` flag per the previous section. 

You can tag your commit for release by running the following bash command:
```bash
git commit -m "Tag for v0.x.x release"  
git tag -a "v0.x.x" -m "Tag message"  
git push origin v0.x.x
```

The final statement will push the tag to Github, where Travis will automatically pick up commit, build and deploy.

> You can check the build status and history on [`travisci/janis`](https://travis-ci.org/PMCC-BioinformaticsCore/janis).
> You can also check [`PyPi`](https://pypi.org/project/janis-pipelines/) to confirm the build has been released.

#### Failing builds for tags

If your tag fails to build, you'll need to:
1. Delete the tag from Github: https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/vX.X.X
2. Delete the tag from your local: `git tag -d 'vX.X.X'`
3. Retag and push per [Tagging and Building](#tagging-and-building)

### Github release notes

You should update the release notes in the `CHANGELOG.md` where you'll need to create a `##` section with three bits of information:

1.  Release commit hash (_eg_: [`bc7f8ad`](https://github.com/PMCC-BioinformaticsCore/janis/commit/bc7f8ad6635c20a8f28788657e8963c2d9ff198e`))
2. Github link to commits comparison: `https://github.com/PMCC-BioinformaticsCore/janis/compare/$prev...$new`
	- Where `$prev` and `$new` is your old, and release commit hash respectively.
3. Release notes: a summary of changes since the last release
