# Janis  (Pre-alpha)


![GitHub stars](https://img.shields.io/github/stars/PMCC-BioinformaticsCore/janis.svg?style=social)  [![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)  [![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest)  [![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines)  [![codecov](https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg)](https://codecov.io/gh/PMCC-BioinformaticsCore/janis) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)
  
  
_Janis is a framework creating specialised, simple workflow definitions that are then transpiled to   
Common Workflow Language or Workflow Definition Language._  
  
Documentation is hosted here: https://janis.readthedocs.io/  
  
## Introduction  

>| WARNING: this project is work-in-progress and is provided as-is without warranty of any kind. There may be breaking changes committed to this repository without notice. |
>|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|


Janis gives you an API to build computational workflows and will generate
a workflow description in CWL and WDL. By using Janis, you get type-safety,
portability and reproducibility across all of your execution environments.


Janis requires a Python installation > 3.6, and can be installed through PIP 
([project page](https://pypi.org/project/janis-pipelines/)):  
  
```bash
# Install janis and the bioinformatics tools
pip3 install janis-pipelines[bioinformatics]  
```  
  
You can import Janis into your project with:  
```python  
import janis as j  
```

## Structure

Janis is structured into components:

- [`core`](https://github.com/PMCC-BioinformaticsCore/janis-core) - contains classes and functions to assist in the workflow building. 
    - [`bioinformatics`](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics) - has bioinformatics tools and data types (`import janis.bioinformatics
`)
    - [`unix`](https://github.com/PMCC-BioinformaticsCore/janis-unix) - Installed by default (`import janis.unix`)
- [`runner`](https://github.com/PMCC-BioinformaticsCore/janis-runner) - An assistant to run Janis workflows using different workflow engines with common semantics.



This repository manages the dependencies for installation and drives the documentation.

### Example  
  
_Further information_: [Simple Workflow](https://janis.readthedocs.io/en/latest/tutorials/echo.html)  
  
Below we've constructed a simple example that takes a string input, calls the 
[echo](https://janis.readthedocs.io/en/latest/tools/unix/echo.html) tool and exposes the 
Echo tool's output as a workflow output.  
  
```python  
import janis as j
from janis.unix.tools.echo import Echo

w = j.Workflow("workflowId")

inp = j.Input("inputIdentifier", data_type=j.String(), value="Hello, World!")
echo = j.Step("stepIdentifier", tool=Echo())
out = j.Output("outputIdentifier")

w.add_edges([
    (inp, echo.inp),  # Connect 'inp' to 'echostep'
    (echo.out, out),  # Connect output of 'echostep' to 'out'
])

# Will print the CWL, input file and relevant tools to the console
w.translate("cwl", to_disk=False)  # or "wdl"

```

We can export a CWL representation to the console using `.translate("cwl")`. By including the 
`to_disk=True` parameter, we can write this workflow to disk at the current location. 
  
#### More examples  

- Bioinformatics workflow tutorial: [AlignSortedBam](https://janis.readthedocs.io/en/latest/tutorials/alignsortedbam.html)
- Unix Toolset: in [`janis/examples`](https://github.com/PMCC-BioinformaticsCore/janis/tree/master/janis/examples).   

- Whole genome germline pipeline: [janis-examplepipelines repository](https://github.com/PMCC-BioinformaticsCore/janis-examplepipelines).  

## Usage

Janis has an API that mirrors the workflow concepts:

- `j.Workflow`: A workflow represents the `Edge`s between `Input`, `Step`, `Output`
  - `j.Input`: An input to a Workflow, has an identifier, a type and a value.
  - `j.Step`: A step also has an identifier and a `Tool` (`CommandTool` or a nested `Workflow`).
  - `j.Output`: An output to a workflow has an identifier and is connected to a step.
  
- `j.CommandTool`: A command line style tool that builds it's command through the inputs and arguments. 
  - `j.ToolInput`: An input to a tool, has an identifier, a type and command line options like `position`, `prefix`
  - `j.ToolArgument`: An argument to a tool that cannot be overridden. Has a value and command line options 
        like `position` and  `prefix`. The value can be a derived type, like an `InputSelector` or `StringFormatter`. 
  - `j.ToolOutput`: Output to a tool, has an identifier, a type and a glob.
  
## About  
  
> _Further information_: [About](https://janis.readthedocs.io/en/latest/about.html)   
  
This project was produced as part of the Portable Pipelines Project in partnership with:    
- [Melbourne Bioinformatics (University of Melbourne) ](https://www.melbournebioinformatics.org.au/)    
- [Peter MacCallum Cancer Centre](https://www.petermac.org/)    
- [Walter and Eliza Hall Institute of Medical Research (WEHI) ](https://www.wehi.edu.au/)    

### References:

Through conference or talks, this project has been referenced by the following titles:

- Walter and Eliza Hall Institute Talk (WEHI) 2019: _Portable Pipelines Project: Developing reproducible bioinformatics pipelines with standardised workflow languages_
- Bioinformatics Open Source Conference (BOSC) 2019: _Janis: an open source tool to machine generate type-safe CWL and WDL workflows_
- Victorian Cancer Bioinformatics Symposium (VCBS) 2019: _Developing portable variant calling pipelines with Janis_
  
  
## Support  
  
### Contributions  
  
> _Further information_: [Development](https://janis.readthedocs.io/en/latest/development/)  
  
This project is work-in-progress and is still in developments. Although we welcome contributions,  due to the immature state of this project we recommend raising issues through the [Github issues page](https://github.com/PMCC-BioinformaticsCore/janis/issues) for Pipeline related issues.  

Information about the project structure and more on contributing can be found within [the documentation](https://janis.readthedocs.io/en/latest/development/).
