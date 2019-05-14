About
######

This project was produced as part of the Portable Pipelines Project in partnership with:

- `Melbourne Bioinformatics (University of Melbourne) <https://www.melbournebioinformatics.org.au/>`_
- `Peter MacCallum Cancer Centre <https://www.petermac.org/>`_
- `Walter and Eliza Hall Institute of Medical Research (WEHI) <https://www.wehi.edu.au//>`_
  

Motivations
===========
  
Given the `awesome list of <https://github.com/pditommaso/awesome-pipeline/>`_ pipeline frameworks, languages and engines, why create another framework to generate workflow langauges.
  
That's a great question, and it's a little complicated. Our project goals are to have a portable workflow specification, that is reproducible across many different compute platforms. And instead of backing one technology, we thought it  would be more powerful to create a technology that can utilise the community's work.  
  
Some additional benefits we get by writing a generic framework is we sanity check connections and also add types that  exist within certain domains. For example within the bioinformatics tools, there's a `BamBai` type that represents an  indexed `.bam` (+ `.bai`) file. With this framework, we don't need to worry about pesky secondary files, or the complications that come when passing them around in WDL either, this framework can take care of that.  

  
  
Related project links:
======================

============== ==================== ====================== ==================== ===============
\              build                docs                   pypi                 codecov
============== ==================== ====================== ==================== ===============
`Janis`_       |Build Status janis| |Documentation Status| |PyPI version janis| |codecov janis|
Bioinformatics |Build Status bio|   See Janis              |PyPI version bio|   \
Shepherd       \                    \                      \                    \
CWL-Gen        |Build Status cwl|   \                      |PyPI version cwl|   |codecov cwl|
WDL-Gen        |Build Status wdl|   \                      |PyPI version wdl|   \
============== ==================== ====================== ==================== ===============

.. _Janis: https://github.com/PMCC-BioinformaticsCore/janis
.. _JanisPIP: https://pypi.org/project/janis-pipelines/

.. |Build Status janis| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis
.. |Documentation Status| image:: https://readthedocs.org/projects/janis/badge/?version=latest
   :target: https://janis.readthedocs.io/en/latest/?badge=latest
.. |PyPI version janis| image:: https://badge.fury.io/py/janis-pipelines.svg
   :target: https://badge.fury.io/py/janis-pipelines
.. |codecov janis| image:: https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/PMCC-BioinformaticsCore/janis
.. |Build Status bio| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics
.. |PyPI version bio| image:: https://badge.fury.io/py/janis-pipelines.bioinformatics.svg
   :target: https://badge.fury.io/py/janis-pipelines.bioinformatics
.. |Build Status cwl| image:: https://travis-ci.org/illusional/python-cwlgen.svg?branch=master
   :target: https://travis-ci.org/common-workflow-language/python-cwlgen
.. |PyPI version cwl| image:: https://badge.fury.io/py/illusional.cwlgen.svg
   :target: https://badge.fury.io/py/illusional.cwlgen
.. |codecov cwl| image:: https://codecov.io/gh/illusional/python-cwlgen/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/illusional/python-cwlgen
.. |Build Status wdl| image:: https://travis-ci.org/illusional/python-wdlgen.svg?branch=master
   :target: https://travis-ci.org/illusional/python-wdlgen
.. |PyPI version wdl| image:: https://badge.fury.io/py/illusional.wdlgen.svg
   :target: https://badge.fury.io/py/illusional.wdlgen
