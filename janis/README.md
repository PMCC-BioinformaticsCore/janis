# Janis
_Portable pipelines assistant_

## Quickstart

Install through PIP ([project page](https://pypi.org/project/janis/)):
```
pip install janis
```

OR

Clone the [GitHub repository](https://github.com/PMCC-BioinformaticsCore/pipelines):
```bash
git clone git@github.com:PMCC-BioinformaticsCore/pipelines.git
```


## About

This project was produced as part of the Portable Pipelines Project in partnership with:
- [Melbourne Bioinformatics (University of Melbourne) ](https://www.melbournebioinformatics.org.au/)
- [Peter MacCallum Cancer Centre](https://www.petermac.org/)
- [Walter and Eliza Hall Institute of Medical Research (WEHI) ](https://www.wehi.edu.au/)

## Usage

You must import `Pipeline` into your project, that is:
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
[Github issues page](https://github.com/PMCC-BioinformaticsCore/pipelines/issues).

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
python -m twine upload dist/portable_pipelines-$VERSION-py3-none-any.whl
``` 

And that's it!
