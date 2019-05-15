

# Janis  

>| WARNING: this project is work-in-progress and is provided as-is without warranty of any kind. There may be breaking changes committed to this repository without notice. |
>|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|

![GitHub stars](https://img.shields.io/github/stars/PMCC-BioinformaticsCore/janis.svg?style=social)  [![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)  [![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest)  [![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines)  [![codecov](https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg)](https://codecov.io/gh/PMCC-BioinformaticsCore/janis)  
  
_Janis is a framework creating specialised, simple workflow definitions that are then transpiled to   
Common Workflow Language or Workflow Definition Language._  
  
Documentation is hosted here: https://janis.readthedocs.io/  
  
## Introduction  
  
Janis is designed to assist in building computational workflows to generate a runnable workflow description (CWL | WDL).  

Janis requires a Python installation > 3.6, and can be installed through PIP ([project page](https://pypi.org/project/janis-pipelines/)) by running:  
  
```bash
# Install janis and the bioinformatics tools
pip3 install janis-pipelines[bioinformatics]  
```  
  
You can import Janis into your project with:  
```python  
import janis as j  
```  
    
### Included tool definitions and types  
  
#### Bioinformatics  
  
The Janis framework can be extended to include a suite of [Bioinformatics data types and tools](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics). These can be installed with the `bioinformatics` install extra option.   
  
```bash  
pip3 install janis-pipelines[bioinformatics]  
```  

#### Unix

Some basic unix tools have been wrapped and included as part of the base Janis module and are the basis for the examples. You can reference these unix tools through `janis.unix.tools`.  
  
These can be referenced by `janis.bioinformatics` or `janis_bioinformatics`, the latter might be easier due to the way  nested python imports work.  
  
## Example  
  
_Further information_: [Simple Workflow](https://janis.readthedocs.io/en/latest/tutorials/simple.html)  
  
Below we've constructed a simple example that takes a string input, uses the [echo](https://janis.readthedocs.io/en/latest/tools/unix/echo.html)   tool to log this to `stdout`, and capturing the `stdout` to output.  to give you a basic idea of how to construct a pipeline.  
  
```python  
import janis as j  
from janis.unix.tools.echo import Echo   
w = j.Workflow("workflow_identifier")  
  
inp = j.Input("input_identifier", j.String())  
step = j.Step("step_identifier", Echo())  
outp = j.Output("output_identifier")  
  
w.add_edges([(inp, step), (step, outp)])
  
# Will print the CWL, input file and relevant tools to the console  
w.translate("cwl")  
```  
We can export a CWL representation to the console using `.translate("cwl")`.   
  
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
  
  
## Support  
  
### Contributions  
  
> _Further information_: [Development](https://janis.readthedocs.io/en/latest/development/)  
  
This project is work-in-progress and is still in developments. Although we welcome contributions,  due to the immature state of this project we recommend raising issues through the [Github issues page](https://github.com/PMCC-BioinformaticsCore/janis/issues) for Pipeline related issues.  
  
If you find an issue with the tool definitions, please see the relevant issue page:  
- [Pipeline-bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues)  
  
Information about the project structure and more on contributing can be found within [the documentation](https://janis.readthedocs.io/en/latest/development/).