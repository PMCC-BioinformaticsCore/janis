# Janis

[![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)
[![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines)

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

Related project links:
- 

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
w.dump_cwl(to_disk=False)
```


## Bioinformatics tools and data types

_Coming soon_

A repository of bioinformatic tools will be build to use within this pipeline. 
The git submodule is embedded here for reference, but can also be found here: [here](https://github.com/PMCC-BioinformaticsCore/pipelines-bioinformatics).

### Intended usage

```
pip install portable-pipeline-bioinformatics
```


## Contributions

Contributions are the bread and butter of open source, and we welcome contributions. 
All sections of this module are written in Python, however a fair understanding of Workflows, CWL or WDL 
might be required to make changes.

If you find an issue with Pipeline related functionality, please report it through the 
[Github issues page](https://github.com/PMCC-BioinformaticsCore/janis/issues).

If you find an issue with the tool definitions, please see the relevant issue page:
- [Pipeline-bioinformatics](https://github.com/PMCC-BioinformaticsCore/pipelines-bioinformatics/issues)


### Releasing Portable Pipelines

Currently the release process is manual. The intent is to run this project through a continuous integration 
system for automatic releases on passing commits to master.

To release, increment the version in `setup.py`, run the following command to generate the release binary
```
python setup.py bdist_wheel
```

To complete the upload, you will need to install `twine` and have your Pip settings configured:

```
python -m twine upload dist/janis_pipelines-$VERSION-py3-none-any.whl
``` 

And that's it!


__________



# WEHI pipeline definition language
A framework for creating specialised, simple workflow definitions that are then converted to Common Workflow Language definitions.

## Python dependencies:

Python dependencies are listed in requirements.txt

Local environment can be initialized by the following command:

>pip install -r requirements.txt

To update the requirements.txt to reflect the latest dependency requirements,

>pip freeze > requirements.txt


Uses NetworkX library: https://networkx.github.io/

## Requirements

This document should consume the WEHI Pipeline Definition Language and emit a translation into CWL (later: WDL). It should contain:
    - A fully specified workflow.
    - A zip of the references (the tools).
    - An input yaml file.

It should have the option to produce a visual graph (as it build the DAG) and perform basic typechecking between connections.

## WEHI Pipeline Definition Language
Following discussions with Evan, I've put together a little _guide_ on how I've interpreted the pipeline language (_name?_).

Theoretically we should be able to consume any input format that will serialize as a dictionary (YAML, JSON), but we'll stick to YAML for the descriptions.

**NB:** There is no more type inferencing if you're seeing this message.

### Linking inputs
You must specify and exactly link inputs to the pipelines, and will likely reference a step's output in another input, you should know what files your tool exports.


```yaml
inputs:                                 # Dictionary of input types, input_labels must be unique
    $input_label:
        $input_type:
            property1: value
            property2: value

    $input_label2:
        $input_type2:
            property1: value
            property2: value

outputs:
    # TBA

steps:    # dictionary of steps
    $step_label:
        tool: $tool/version             # Should be able to just say 'tool' or 'toolCategory' as well
        inputs:
            $input1: $input_label
            $input2: $input_label2

    $step_label2:
        tool: $tool_category
        $input1: $step_label/output1    # You must refer to a tool's documentation to find out the types it exports
        $input2: $input_label2
```

For a very simplified step, I'd be happy to just specify the tool, and all outputs from the previous step is provided to the current step, ie. a purely linear workflow, eg:
```yaml
steps:
    $step1: $tool1
    $step2: $tool2
    $step3: $tool3
```

Notes:
- You must NOT have a label called "input"

### Questions

- Inheritance with type checking
- How are secondary files dealt with, are they all passed?
- Should inputs include actual values to generate the yaml?