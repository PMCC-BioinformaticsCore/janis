About
######

This project was produced as part of the Portable Pipelines Project in partnership with:

- `Melbourne Bioinformatics (University of Melbourne) <https://www.melbournebioinformatics.org.au/>`_
- `Peter MacCallum Cancer Centre <https://www.petermac.org/>`_
- `Walter and Eliza Hall Institute of Medical Research (WEHI) <https://www.wehi.edu.au//>`_
  

Motivations
===========
  
Given the `awesome list of <https://github.com/pditommaso/awesome-pipeline/>`_ pipeline frameworks, languages and engines, Janis was created to facilitate better interop and accessibility within the workflow space. 

Janis' translate ability is designed as a productivity tool which assists in migrating tools / workflows from one specification to another. 

The intended use-cases include:

- Research groups migrating their pipelines to a new specification
- Migrating a published workflow to a specification you are more familiar with
- Migrating a workflow to another specification for execution on a specific compute infrastructure
  

Components
==========

Janis consists of multiple modules which can be used independently or together.

| - **Janis Core**
|   The core python API and translate functionality of Janis. Contains the core datatypes, tool, and workflow model definitions. 
|
| - **Janis Bioinformatics**
|   Datatypes, as well as implemented Tools & Workflows for bioinformatics pipelines.
| 
| - **Janis Unix**
|   Datatypes, as well as implemented Tools & Workflows for generic unix software pipelines. 
|
| - **Janis Assistant**
|   As part of the Portable Pipelines Project, we produced an execution assistant called ``janis-assistant``. 
|   Its purpose is to run workflows written in Janis, track the progress and report the results back in a robust way.



  
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
