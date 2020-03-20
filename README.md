# Janis (janis-pipelines) (Alpha)


[![GitHub stars](https://img.shields.io/github/stars/PMCC-BioinformaticsCore/janis.svg?style=social)](https://github.com/PMCC-BioinformaticsCore/janis) [![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)  [![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest)  [![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines)  [![codecov](https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg)](https://codecov.io/gh/PMCC-BioinformaticsCore/janis) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black) [![Gitter chat](https://badges.gitter.im/janis-pipelines.png)](https://gitter.im/janis-pipelines/community)
  
_Janis is a framework creating specialised, simple workflow definitions that are then transpiled to   
Common Workflow Language or Workflow Definition Language._  
  
Documentation is available here: https://janis.readthedocs.io/

## Introduction  

>| WARNING: this project is work-in-progress and is provided as-is without warranty of any kind. There may be breaking changes committed to this repository without notice. |
>|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------|


Janis gives you an API to build computational workflows and will generate
a workflow description in CWL and WDL. By using Janis, you get type-safety,
portability and reproducibility across all of your execution environments.


Janis requires a Python installation > 3.6, and can be installed through PIP 
([project page](https://pypi.org/project/janis-pipelines/)):  
  
```bash
# Install janis and the toolkits
pip3 install janis-pipelines 
```

There are two ways to use Janis:

- Build workflows (and translate to CWL / WDL)
- Run tools or workflows with CWLTool or Cromwell

### Example workflow

  
Let's construct a simple example that takes a string input, calls the 
[echo](https://janis.readthedocs.io/en/latest/tools/unix/echo.html) tool and exposes the 
Echo tool's output as a workflow output. 


  
```bash
# write the workflow to `helloworld.py`
cat <<EOT >> helloworld.py
import janis as j
from janis_unix.tools import Echo

w = j.WorkflowBuilder("hello_world")

w.input("input_to_print", j.String)
w.step("echo", Echo(inp=w.input_to_print))
w.output("echo_out", source=w.echo.out)
EOT


# Translate workflow to WDL
janis translate helloworld.py wdl

# Run the workflow
janis run -o helloworld-tutorial helloworld.py --input_to_print "Hello, World!"

# See your output
cat helloworld-tutorial/echo_out
# Hello, World!
```

  
### How to use Janis

- [Tutorial 0 - Introduction to Janis](https://janis.readthedocs.io/en/latest/tutorials/tutorial0.html)
- [Tutorial 1 - Building a workflow](https://janis.readthedocs.io/en/latest/tutorials/tutorial1.html)
- [Tutorial 2 - Wrapping a new tool](https://janis.readthedocs.io/en/latest/tutorials/tutorial2.html)


#### Workshops

In addition, there are fully self-guided workshops that more broadly go through the functionality of Janis:

- [Workshop 1](https://github.com/PMCC-BioinformaticsCore/janis-workshops/tree/master/workshop1)
- [Workshop 2](https://github.com/PMCC-BioinformaticsCore/janis-workshops/tree/master/workshop2)

### Examples

Sometimes it's easier to learn by examples, here are a few hand picked examples:

- [Samtools View](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/blob/master/janis_bioinformatics/tools/samtools/view/base.py) ([Docs](https://janis.readthedocs.io/en/latest/tools/bioinformatics/samtools/samtoolsview.html))

- [WGS Germline pipeline (GATK Only)](https://github.com/PMCC-BioinformaticsCore/janis-pipelines/blob/master/janis_pipelines/wgs_germline_gatk/wgsgermlinegatk.py) ([Docs](https://janis.readthedocs.io/en/latest/pipelines/wgsgermlinegatk.html))


### Toolbox

There are two toolboxes currently available on Janis:

- [Unix](https://github.com/PMCC-BioinformaticsCore/janis-unix) ([list of tools](https://janis.readthedocs.io/en/latest/tools/bioinformatics/index.html))
- [Bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics) ([list of tools](https://janis.readthedocs.io/en/latest/tools/unix/index.html))



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
- GIW / ABACBS 2019: _Janis: A Python framework for Portable Pipelines_
- Australian BioCommons, December 2019: _Portable pipelines: build once and run everywhere with Janis_
  
  
## Support  

### v0.9.0 Backwards Compatability

**NOTE: Version 0.9.0 brings changes to output directories and camel case changes**

- Janis watch will be incompatible with previously run workflows
- Your configs might break, as previous versions of janis were not cautious about camel case.
- Your templates might not work with unrecognised keys (try changing them to camel case instead)
- Changes to BamBai indexes, format is now `.bam.bai`

See the [CHANGELOG](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/CHANGELOG.md)
for more information.


### Contributions
  
> _Further information_: [Development](https://janis.readthedocs.io/en/latest/development/)  
  
To get help with Janis, please ask a question on [Gitter](ttps://gitter.im/janis-pipelines/community) or 
[raise an issue](https://github.com/PMCC-BioinformaticsCore/janis/issues) on GitHub.

If you find an issue with the tool definitions, please see the relevant issue page:

- [Pipeline-bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues)

This project is work-in-progress and is still in developments. Although we welcome contributions,
due to the immature state of this project we recommend raising issues through the
[Github issues page](https://github.com/PMCC-BioinformaticsCore/janis/issues) for Pipeline related issues.

Information about the project structure and more on contributing can be found within the documentation.
