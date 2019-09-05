Welcome to Janis!
********************************

.. image:: https://img.shields.io/github/stars/PMCC-BioinformaticsCore/janis.svg?style=social
    :target: https://github.com/PMCC-BioinformaticsCore/janis
    :alt: GitHub stars

.. image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master&style=flat
    :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis
    :alt: Travis Build Status

.. image:: https://readthedocs.org/projects/janis/badge/?version=latest
    :target: https://janis.readthedocs.io/en/latest/?badge=latest)
    :alt: Documentation

.. image:: https://badge.fury.io/py/janis-pipelines.svg
    :target: https://pypi.org/project/janis-pipelines/
    :alt: Pypi module

.. image:: https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/PMCC-BioinformaticsCore/janis
    :alt: Code Coverage

.. note::
	This project is *work-in-progress* and is provided as-is without warranty of any kind. There may be breaking changes
	committed to this repository without notice.

Janis is a framework creating specialised, simple workflow definitions that are then transpiled to
Common Workflow Language or Workflow Definition Language.

It was developed as part of the Portable Pipelines Project, a collaboration between:

- `Melbourne Bioinformatics (University of Melbourne) <https://www.melbournebioinformatics.org.au/>`_
- `Peter MacCallum Cancer Centre <https://www.petermac.org/>`_ 
- `Walter and Eliza Hall Institute of Medical Research (WEHI) <https://www.wehi.edu.au/>`_.


Introduction
============

Janis is designed to assist in building computational workflows to generate a runnable workflow description (CWL | WDL).
It can be installed through the Janis `project page <https://pypi.org/project/janis-pipelines/>`__ by running:

.. code-block:: bash

   pip3 install janis-pipelines

You can import Janis into your project by:

.. code-block:: python

   import janis as j

Included tool definitions and types
===================================

Some basic unix tools have been wrapped and included as part of the base Janis module and are the basis for the examples.
You can reference these unix tools through `janis.unix.tools`.

Bioinformatics
--------------

The Janis framework can be extended to include a suite of
`Bioinformatics data types and tools <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`__. These can be
installed with the ``bioinformatics`` install extra option.

.. code-block:: bash

   pip3 install janis-pipelines[bioinformatics]

These can be referenced by ``janis.bioinformatics`` or ``janis_bioinformatics``, the latter might be easier due to the way
nested python imports work.


Example
========

*Further information*: `Simple Workflow </tutorials/simple.html>`__

Below we've constructed a simple example that takes a string input, uses the `echo <tools/unix/echo.html>`__)
tool to log this to ``stdout``, and explicitly outputting this ``stdout`` to give you a basic idea of how to construct a pipeline.

.. code-block:: python

   import janis as j
   from janis.unix.tools.echo import Echo

   w = j.Workflow("workflowId")

   w.input("inputIdentifier", j.String, default="Hello, World!")
   w.step("stepIdentifier", Echo, inp=w.inputIdentifier)
   w.output("outputIdentifier", source=w.stepIdentifier.out)

   # Will print the CWL, input file and relevant tools to the console
   w.translate("cwl", to_disk=False)  # or "wdl"

More information can be found on creating edges on the
`Building Connections <tutorials/buildingconnections.html>`__ documentation.

Now we've created our workflow, we can export a CWL representation to the console using ``.translate("cwl")``.

More examples
-------------

There are some simple example pipelines that use the unix toolset in
`janis/examples <https://github.com/PMCC-BioinformaticsCore/janis/tree/master/janis/examples>`__.

Additionally there are example bioinformatics workflows that use Janis and the bioinformatics tools in the
`janis-examplepipelines repository <https://github.com/PMCC-BioinformaticsCore/janis-pipelines>`__.


Support
=======

This project is work-in-progress and is still in developments. Although we welcome contributions,
due to the immature state of this project we recommend raising issues through the
`Github issues page <https://github.com/PMCC-BioinformaticsCore/janis/issues>`__ for Pipeline related issues.

If you find an issue with the tool definitions, please see the relevant issue page:

* `Pipeline-bioinformatics <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues>`__

Information about the project structure and more on [contributing]() can be found within the documentation.

Releasing Janis
---------------

   *Further Information*: `Releasing <development/releasing.html>`__

Releasing is automatic! Simply increment the version number in `setup.py` (`SemVer <https://semver.org>`__),
and tag that commit with the same version identifier:

.. code-block:: bash
   git commit -m "Tag for v0.x.x release"
   git tag -a "v0.x.x" -m "Tag message"
   git push --follow-tags


Dependencies and Resources
==========================

Janis includes a number of dependencies

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


Contents
========

.. toctree::
   :maxdepth: 1
   :hidden:

   self
   about
   userguide

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Tutorials

   tutorials/gettingstarted
   tutorials/simple
   tutorials/alignsortedbam
   tutorials/running

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: References

   references/index
   tools/index
   datatypes/index


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Development

   development/index
   development/testing
   development/releasing
   development/faq

Indices and tables
^^^^^^^^^^^^^^^^^^

* :ref:`genindex`
* :ref:`search`
.. * :ref:`modindex`
