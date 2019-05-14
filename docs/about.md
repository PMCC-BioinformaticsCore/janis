# About

This project was produced as part of the Portable Pipelines Project in partnership with:  
- [Melbourne Bioinformatics (University of Melbourne) ](https://www.melbournebioinformatics.org.au/)  
- [Peter MacCallum Cancer Centre](https://www.petermac.org/)  
- [Walter and Eliza Hall Institute of Medical Research (WEHI) ](https://www.wehi.edu.au/)  
  
  
### Motivations  
  
Given the [awesome list of](https://github.com/pditommaso/awesome-pipeline) pipeline frameworks, languages and engines, why create another framework to generate workflow langauges.  
  
That's a great question, and it's a little complicated. Our project goals are to have a portable workflow specification, that is reproducible across many different compute platforms. And instead of backing one technology, we thought it  would be more powerful to create a technology that can utilise the community's work.  
  
Some additional benefits we get by writing a generic framework is we sanity check connections and also add types that  exist within certain domains. For example within the bioinformatics tools, there's a `BamBai` type that represents an  indexed `.bam` (+ `.bai`) file. With this framework, we don't need to worry about pesky secondary files, or the complications that come when passing them around in WDL either, this framework can take care of that.  

 
  
  
### Related project links:  

  
|  | build  | docs  | pypi | codecov |  
|---|:-:|:-:|:-:|:-:|  
| [Janis]([https://github.com/PMCC-BioinformaticsCore/janis](https://github.com/PMCC-BioinformaticsCore/janis)) |  [![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis)  | [![Documentation Status](https://readthedocs.org/projects/janis/badge/?version=latest)](https://janis.readthedocs.io/en/latest/?badge=latest) | [![PyPI version](https://badge.fury.io/py/janis-pipelines.svg)](https://badge.fury.io/py/janis-pipelines) |  [![codecov](https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg)](https://codecov.io/gh/PMCC-BioinformaticsCore/janis) |  
| [Bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics) |[![Build Status](https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics.svg?branch=master)](https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics) | See Janis  |  [![PyPI version](https://badge.fury.io/py/janis-pipelines.bioinformatics.svg)](https://badge.fury.io/py/janis-pipelines.bioinformatics)  |   |  
| [Shepherd](https://github.com/PMCC-BioinformaticsCore/shepherd) | | | | |  
| [CWL-Gen (fork)](https://github.com/illusional/python-cwlgen) | [![Build Status](https://travis-ci.org/illusional/python-cwlgen.svg?branch=master)](https://travis-ci.org/illusional/python-cwlgen) |   |[![PyPI version](https://badge.fury.io/py/illusional.cwlgen.svg)](https://badge.fury.io/py/illusional.cwlgen) | [![codecov](https://codecov.io/gh/illusional/python-cwlgen/branch/master/graph/badge.svg)](https://codecov.io/gh/illusional/python-cwlgen)|  
| [WDL-Gen](https://github.com/illusional/python-wdlgen) | [![Build Status](https://travis-ci.org/illusional/python-wdlgen.svg?branch=master)](https://travis-ci.org/illusional/python-wdlgen) || [![PyPI version](https://badge.fury.io/py/illusional.wdlgen.svg)](https://badge.fury.io/py/illusional.wdlgen) | |