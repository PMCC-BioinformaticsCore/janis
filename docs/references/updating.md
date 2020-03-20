# Updating Janis

Janis may seem like (an _organised_) chaos of individual Python packages, but they're individual to ensure proper separation of concerns, and allows us to provide bug fixes to each section separately.

The primary Janis module `janis-pipelines` ([RELEASES](https://github.com/PMCC-BioinformaticsCore/janis/releases) page) contains a list of all the packages and their versions. But sometimes we make quick bug fixes to individual projects, without doing a big release of the whole project.

We recommend sticking to these major releases, you can update Janis to the latest major release with:

```bash
pip3 install --no-cache --upgrade janis-pipelines
```

## Power users

Okay, so you're a power user and like to live on the edge. You need to decide now whether you like to live close to the cliff, or on the bleeding edge:

- Close-ish: Installing specific modules from Pip
- Bleeding edge: Installing individual projects from GitHub.

### Installing from Pip

Janis projects follow a predictable naming structure:

> With one except, janis-assistant is called `janis-pipelines.runner` for legacy reasons.

```
janis-pipelines.<projectname>
```

You could for example update janis-core with:

```
pip3 install --no-cache --upgrade janis-pipelines.core
```

### Bleeding edge: From GitHub

You should be warned that builds from GitHub are not always stable, and sometimes projects get released together as they contain inter-dependencies.
We'd also recommend that you install from GitHub with the `--no-dependencies` flag.

Due diligence over, you can install from GitHub with the following line (replacing the GitHub link with the GitHub repo you're trying to install).

```bash
pip3 install --no-cache --upgrade --no-dependencies git+https://github.com/PMCC-BioinformaticsCore/janis-<project>.git
```


