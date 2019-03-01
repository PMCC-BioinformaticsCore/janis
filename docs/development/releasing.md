# Releasing the Portable Pipeline Projects


_This page is under construction_

## Janis:

Releasing is automatic! Simply increment the version number in `setup.py` ([SemVer](https://semver.org)), 
and tag that commit with the same version identifier:
```bash
git commit -m "Tag for v0.x.x release"
git tag -a "v0.x.x" -m "Tag message"
git push --follow-tags
```

[Travis](https://travis-ci.org/PMCC-BioinformaticsCore/janis) will automatically build the repository, 
run the associated unit tests and if they succeed, deploy to [PyPi](https://pypi.org/project/janis-pipelines/). 