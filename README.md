# Janis

> | WARNING: this project is work-in-progress and is provided as-is without warranty of any kind. There may be breaking changes committed to this repository without notice.  |
> | :-- |

![GitHub stars](https://img.shields.io/github/stars/PMCC-BioinformaticsCore/janis.svg?style=social)
[![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)
[![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines)
[![codecov](https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg)](https://codecov.io/gh/PMCC-BioinformaticsCore/janis)

_Janis is a framework creating specialised, simple workflow definitions that are then transpiled to 
Common Workflow Language or Workflow Definition Language._

Documentation is hosted here: https://janis.readthedocs.io/

## Introduction

Janis is designed to assist in building computational workflows to generate a runnable workflow description (CWL | WDL).
It can be installed through PIP ([project page](https://pypi.org/project/janis-pipelines/)) by running:

```
pip install janis-pipelines
```

You can import Janis into your project by:
```python
import janis as j
```

### Included tool definitions and types

Some basic unix tools have been wrapped and included as part of the base Janis module and are the basis for the examples.
You can reference these unix tools through `janis.unix.tools`.

#### Bioinformatics

The Janis framework can be extended to include a suite of 
[Bioinformatics data types and tools](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics). These can be
installed with the `bioinformatics` install extra option. 

`pip install janis-pipelines[bioinformatics]`

These can be referenced by `janis.bioinformatics` or `janis_bioinformatics`, the latter might be easier due to the way
nested python imports work.

### Example

_Further information_: [Simple Workflow](https://janis.readthedocs.io/en/latest/tutorials/simple.html)

Below we've constructed a simple example that takes a string input, uses the [echo](https://janis.readthedocs.io/en/latest/tools/unix/echo.html) 
tool to log this to `stdout`, and explicitly outputting this `stdout` to give you a basic idea of how to construct a pipeline.

```python
import janis as j
from janis.unix.tools.echo import Echo 

w = j.Workflow("workflow_identifier")

inp = j.Input("input_identifier", j.String())
step = j.Step("step_identifier", Echo())
outp = j.Output("output_identifier")

w.add_pipe(inp, step, outp)

# Will print the CWL, input file and relevant tools to the console
w.dump_translation("cwl")
```

The `add_pipe` method is aware of the inputs and outputs of the arguments you provide it, and automatically
joins the relevant non-optional parts together. More information can be found on creating edges on the 
["Building Connections"](https://janis.readthedocs.io/en/latest/tutorials/buildingconnections.html) documentation.

Now we can an in-memory workflow, we can export a CWL representation to the console using `.dump_translation("cwl")`. 

#### More examples

There are some simple example pipelines that use the unix toolset in 
[`janis/examples`](https://github.com/PMCC-BioinformaticsCore/janis/tree/master/janis/examples). 

Additionally there are example bioinformatics workflows that use Janis and the bioinformatics tools in the 
[janis-examplepipelines repository](https://github.com/PMCC-BioinformaticsCore/janis-examplepipelines).


## About

This project was produced as part of the Portable Pipelines Project in partnership with:
- [Melbourne Bioinformatics (University of Melbourne) ](https://www.melbournebioinformatics.org.au/)
- [Peter MacCallum Cancer Centre](https://www.petermac.org/)
- [Walter and Eliza Hall Institute of Medical Research (WEHI) ](https://www.wehi.edu.au/)


### Motivations

Given the [awesome list of](https://github.com/pditommaso/awesome-pipeline) pipeline frameworks, languages and engines,
why create another framework to generate workflow langauges.

That's a great question, and it's a little complicated. Our project goals are to have a portable workflow specification,
that is reproducible across many different compute platforms. And instead of backing one technology, we thought it 
would be more powerful to create a technology that can utilise the community's work.

Some additional benefits we get by writing a generic framework is we sanity check connections and also add types that 
exist within certain domains. For example within the bioinformatics tools, there's a `BamBai` type that represents an 
indexed `.bam` (+ `.bai`) file. With this framework, we don't need to worry about pesky secondary files, or the
complications that come when passing them around in WDL either, this framework can take care of that.


### Related project links:
- Janis:
    - Janis Documentation: https://janis.readthedocs.io/en/latest/
    - Janis Git: https://github.com/PMCC-BioinformaticsCore/janis
    - Janis PyPi: https://pypi.org/project/janis-pipelines/
    - Janis Bioinformatics: 
        - Git: https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics
        - PyPi: https://pypi.org/project/illusional.wdlgen/

- Shepherd: https://github.com/PMCC-BioinformaticsCore/shepherd

- CWLGen (forked): https://github.com/illusional/python-cwlgen
- WDLGen: https://github.com/illusional/python-wdlgen

|  | build  | docs  | pypi | codecov |
|---|:-:|:-:|:-:|:-:|
| Janis |  [![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)  | [![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest) | [![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines) |  [![codecov](https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg)](https://codecov.io/gh/PMCC-BioinformaticsCore/janis) |
| Bioinformatics |[![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics) | See Janis  |  [![PyPI version](https://badge.fury.io/py/janis-pipelines.bioinformatics.svg)](https://badge.fury.io/py/janis-pipelines.bioinformatics)  |   |
| Shepherd | | | | |
| CWL-Gen | [![Build Status](https://travis-ci.org/illusional/python-cwlgen.svg?branch=master)](https://travis-ci.org/common-workflow-language/python-cwlgen) |   |[![PyPI version](https://badge.fury.io/py/illusional.cwlgen.svg)](https://badge.fury.io/py/illusional.cwlgen) | [![codecov](https://codecov.io/gh/illusional/python-cwlgen/branch/master/graph/badge.svg)](https://codecov.io/gh/illusional/python-cwlgen)|
| WDL-Gen | [![Build Status](https://travis-ci.org/illusional/python-wdlgen.svg?branch=master)](https://travis-ci.org/illusional/python-wdlgen) || [![PyPI version](https://badge.fury.io/py/illusional.wdlgen.svg)](https://badge.fury.io/py/illusional.wdlgen) | |


## Support

### Contributions

This project is work-in-progress and is still in developments. Although we welcome contributions,
due to the immature state of this project we recommend raising issues through the 
[Github issues page](https://github.com/PMCC-BioinformaticsCore/janis/issues) for Pipeline related issues.

If you find an issue with the tool definitions, please see the relevant issue page:
- [Pipeline-bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues)

Information about the project structure and more on [contributing]() can be found within the documentation.

### Releasing Janis

> _Further Information_: [Releasing](https://janis.readthedocs.io/en/latest/development/releasing.html)

Releasing is automatic! Simply increment the version number in `setup.py` ([SemVer](https://semver.org)), 
and tag that commit with the same version identifier:
```bash
git commit -m "Tag for v0.x.x release"
git tag -a "v0.x.x" -m "Tag message"
git push --follow-tags
```