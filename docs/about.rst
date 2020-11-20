About
######

This project was produced as part of the Portable Pipelines Project in partnership with:

- `Melbourne Bioinformatics (University of Melbourne) <https://www.melbournebioinformatics.org.au/>`_
- `Peter MacCallum Cancer Centre <https://www.petermac.org/>`_
- `Walter and Eliza Hall Institute of Medical Research (WEHI) <https://www.wehi.edu.au//>`_
  

Motivations
===========
  
Given the `awesome list of <https://github.com/pditommaso/awesome-pipeline/>`_ pipeline frameworks, languages and engines, why create another framework to generate workflow languages?
  
That's a great question, and it's a little complicated. Our project goals are to have a portable workflow specification, that is reproducible across many different compute platforms. And instead of backing one technology, we thought it  would be more powerful to create a technology that can utilise the community's work.  
  
Some additional benefits we get by writing a generic framework is we sanity check connections and also add types that  exist within certain domains. For example within the bioinformatics tools, there's a `BamBai` type that represents an  indexed `.bam` (+ `.bai`) file. With this framework, we don't need to worry about pesky secondary files, or the complications that come when passing them around in WDL either, this framework can take care of that.  


Assistant
=========

As part of the Portable Pipelines Project, we produced an execution assistant called ``janis-assistant``. Its purpose is
to run workflows written in Janis, track the progress and report the results back in a robust way.


  
Related project links:
======================

======================= ====================== ====================== ==================== ===============
\                       build                  docs                   pypi                 codecov
======================= ====================== ====================== ==================== ===============
`Janis`_                |Build Status janis|   |Documentation Status| |PyPI version janis| Docs only
`Janis core`_           |Build Status core|    See Janis              |PyPI version core|  |codecov core|
`Janis bioinformatics`_ |Build Status bio|     See Janis              |PyPI version bio|   \
`Janis unix`_           |Build Status unix|    See Janis              |PyPI version unix|  \
`Janis assistant`_      |Build Status asst|    See Janis              |PyPI version asst|  \
======================= ====================== ====================== ==================== ===============

.. _Janis: https://github.com/PMCC-BioinformaticsCore/janis
.. _Janis core: https://github.com/PMCC-BioinformaticsCore/janis-core
.. _Janis bioinformatics: https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics
.. _Janis unix: https://github.com/PMCC-BioinformaticsCore/janis-unix
.. _Janis assistant: https://github.com/PMCC-BioinformaticsCore/janis-assistant

.. _JanisPIP: https://pypi.org/project/janis-pipelines/

.. |Documentation Status| image:: https://readthedocs.org/projects/janis/badge/?version=latest
   :target: https://janis.readthedocs.io/en/latest/?badge=latest

.. |Build Status janis| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis
.. |Build Status core| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis-core.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis-core
.. |Build Status bio| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis-bioinformatics
.. |Build Status unix| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis-unix.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis-unix
.. |Build Status asst| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis-runner.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis-runner


.. |PyPI version janis| image:: https://badge.fury.io/py/janis-pipelines.svg
   :target: https://badge.fury.io/py/janis-pipelines
.. |PyPI version core| image:: https://badge.fury.io/py/janis-pipelines.core.svg
   :target: https://badge.fury.io/py/janis-pipelines.core
.. |PyPI version unix| image:: https://badge.fury.io/py/janis-pipelines.unix.svg
   :target: https://badge.fury.io/py/janis-pipelines.unix
.. |PyPI version bio| image:: https://badge.fury.io/py/janis-pipelines.bioinformatics.svg
   :target: https://badge.fury.io/py/janis-pipelines.bioinformatics
.. |PyPI version asst| image:: https://badge.fury.io/py/janis-pipelines.runner.svg
   :target: https://badge.fury.io/py/janis-pipelines.runner



.. |codecov core| image:: https://codecov.io/gh/PMCC-BioinformaticsCore/janis-core/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/PMCC-BioinformaticsCore/janis-core
