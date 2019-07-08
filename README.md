

# Janis  


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

### Example  
  
_Further information_: [Simple Workflow](https://janis.readthedocs.io/en/latest/tutorials/echo.html)  
  
Below we've constructed a simple example that takes a string input, calls the 
[echo](https://janis.readthedocs.io/en/latest/tools/unix/echo.html) tool and exposes the 
Echo tool's output as a workflow output.  
  
```python  
import janis as j  
from janis.unix.tools.echo import Echo   

w = j.Workflow("workflowId")  
  
inp = j.Input("inputIdentifier", j.String(), value="my value to print")  
echostep = j.Step("stepIdentifier", Echo())  
outp = j.Output("outputIdentifier")  
  

w.add_edges([
    (inp, echostep.inp),    # Connect 'inp' to 'echostep'
    (echostep, outp.outp)    # Connect output of 'echostep' to 'out'
])
  
# Will print the CWL, input file and relevant tools to the console  
w.translate("cwl")  # or "wdl"
```

We can export a CWL representation to the console using `.translate("cwl")`.

#### Named inputs and Outputs

Every input and output of a tool is named. In this example, Janis knows that there is only one
input and one output of `echostep`, so can automatically connect these together. You should see
a statement in the console that indicates that Janis has automatically made this connection.

```
[INFO]: The node 'stepIdentifier' was not a fully qualified input of the tool 'Echo', this was automatically corrected (stepIdentifier → stepIdentifier.inp)
[INFO]: The node 'outputIdentifier' under-referenced an output the step 'stepIdentifier' (tool: 'Echo'), this was automatically corrected (stepIdentifier → stepIdentifier.outp)
```

    
### Included tool definitions and types  
  
#### Bioinformatics  
  
The Janis framework can be extended to include a suite of 
[Bioinformatics data types and tools](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics). 
These can be installed with the `bioinformatics` install extra option.   
  
```bash  
pip3 install janis-pipelines[bioinformatics]  
```  

#### Unix

> Tool document: 

Some basic unix tools have been wrapped and included as part of the base Janis module and 
are the basis for the examples. You can reference these unix tools through 
`janis.unix.tools`.  
  
These can be referenced by `janis.bioinformatics` or `janis_bioinformatics`, the latter might be easier due to the way  nested python imports work.  
  
   
  
#### More examples  

- Bioinformatics workflow tutorial: [AlignSortedBam](https://janis.readthedocs.io/en/latest/tutorials/alignsortedbam.html)
- Unix Toolset: in [`janis/examples`](https://github.com/PMCC-BioinformaticsCore/janis/tree/master/janis/examples).   

- Whole genome germline pipeline: [janis-examplepipelines repository](https://github.com/PMCC-BioinformaticsCore/janis-examplepipelines).  
  
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
  
If you find an issue with the tool definitions, please see the relevant issue page:  
- [Pipeline-bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues)  
  
Information about the project structure and more on contributing can be found within [the documentation](https://janis.readthedocs.io/en/latest/development/).