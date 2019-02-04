# Janis

[![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)
[![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines)
[![codecov](https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg)](https://codecov.io/gh/PMCC-BioinformaticsCore/janis)

_Portable pipelines assistant_: A framework for creating specialised, simple workflow definitions that are then converted to Common Workflow Language or Workflow Definition Language.


## Introduction

Janis is designed to assist in building pipelines. It has a collection of prebuilt tools

Install through PIP ([project page](https://pypi.org/project/janis-pipelines/)):
```
pip install janis-pipelines
```

And ReadTheDocs: https://janis.readthedocs.io/en/latest/
___
OR

Clone the [GitHub repository](https://github.com/PMCC-BioinformaticsCore/janis):
```bash
git clone git@github.com:PMCC-BioinformaticsCore/janis.git
```

## About

This project was produced as part of the Portable Pipelines Project in partnership with:
- [Melbourne Bioinformatics (University of Melbourne) ](https://www.melbournebioinformatics.org.au/)
- [Peter MacCallum Cancer Centre](https://www.petermac.org/)
- [Walter and Eliza Hall Institute of Medical Research (WEHI) ](https://www.wehi.edu.au/)

### Related project links:
- Janis:
    - Janis Documentation: https://janis.readthedocs.io/en/latest/
    - Janis Git: https://github.com/PMCC-BioinformaticsCore/janis
    - Janis PyPi: https://pypi.org/project/janis-pipelines/
    - Janis Bioinformatics: 
        - Git: https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics
        - PyPi: _Coming soon_1

- CWLGen (forked): https://github.com/illusional/python-cwlgen
- WDLGen: https://github.com/illusional/python-wdlgen



## Usage

You must import `janis` into your project, that is:
```python
import janis as j
``` 

### Included definitions

Some unix tools have been wrapped and included as part of the pip module. They are located at `Pipeline.unix.tools/`.
The examples will use the included unix tools, with more information about bioinformatics tools down below. 
See the section about contributions if you find an error in the tool definitions.

### Creating workflows

A Workflow consists of inputs, outputs and steps (which each have their own tool).
You can connect these components together with edges. Let's look the simple untar workflow.

```python
import janis as j
from janis.unix.tools.echo import Echo 

w = j.Workflow("workflow_identifier")

inp = j.Input("input_identifier", p.String())
step = j.Step("step_identifier", Echo())
outp = j.Output("output_identifier")

w.add_pipe(inp, step, outp)

# Will print the CWL, input file and relevant tools to the console
w.dump_translation("cwl", to_disk=False)
```


## Bioinformatics tools and data types

_Coming soon_

A repository of bioinformatic tools will be build to use within this pipeline. 
The git submodule is embedded here for reference, but can also be found here: [here](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics).

### Intended usage

```
pip install janis-pipelines[bioinformatics]
```

Then you can simple import:
```
import janis.bioinformatics
```


## Contributions

Contributions are the bread and butter of open source, and we welcome contributions. 
All sections of this module are written in Python, however a fair understanding of Workflows, CWL or WDL 
might be required to make changes.

If you find an issue with Pipeline related functionality, please report it through the 
[Github issues page](https://github.com/PMCC-BioinformaticsCore/janis/issues).

If you find an issue with the tool definitions, please see the relevant issue page:
- [Pipeline-bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues)


### Releasing Portable Pipelines

Releasing is automatic! Simply increment the version number in `setup.py` ([SemVer](https://semver.org)), 
and tag that commit with the same version identifier:
```
git tag -a "v0.x.x" -m "Tag message"
git push --follow-tags
```

[Travis](https://travis-ci.org/PMCC-BioinformaticsCore/janis) will automatically build the repository, 
run the associated unit tests and if they succeed, deploy to [PyPi](https://pypi.org/project/janis-pipelines/). 
