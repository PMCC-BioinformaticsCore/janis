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

.. image:: https://badges.gitter.im/janis-pipelines.png
    :target: https://gitter.im/janis-pipelines/community
    :alt: Gitter chat

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

A handful of unix and bioinformatics tools have been included in your installation with Janis.

Bioinformatics
--------------

The Janis framework can be extended to include a suite of
`Bioinformatics data types and tools <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`__. These are installed by default.

This can be referenced through ``janis_bioinformatics``, the tool documentation suggests how this may be imported.

Example
========

Below we've constructed a simple example that takes a string input, uses the `echo <tools/unix/echo.html>`_
tool to log this to ``stdout``, and explicitly outputting this ``stdout`` to give you a basic idea of how to construct a pipeline.

.. code-block:: python

   import janis as j
   from janis.tools import Echo

   w = j.WorkflowBuilder("workflowId")

   w.input("inputIdentifier", j.String, default="Hello, World!")
   w.step("stepIdentifier", Echo(inp=w.inputIdentifier))
   w.output("outputIdentifier", source=w.stepIdentifier.out)

   # Will print the CWL, input file and relevant tools to the console
   w.translate("cwl", to_disk=False)  # or "wdl"


Now we've created our workflow, we can export a CWL representation to the console using ``.translate("cwl")``.

More examples
-------------

There are some simple example pipelines that use the unix toolset in
`janis/examples <https://github.com/PMCC-BioinformaticsCore/janis/tree/master/janis/examples>`__.

Additionally there are example bioinformatics workflows that use Janis and the bioinformatics tools in the
`janis-pipelines repository <https://github.com/PMCC-BioinformaticsCore/janis-pipelines>`__.


Support
=======

To get help with Janis, please ask a question on `Gitter <https://gitter.im/janis-pipelines/community>`__ or 
`raise an issue <https://github.com/PMCC-BioinformaticsCore/janis/issues>`__ on GitHub.

If you find an issue with the tool definitions, please see the relevant issue page:

* `Pipeline-bioinformatics <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues>`__

This project is work-in-progress and is still in developments. Although we welcome contributions,
due to the immature state of this project we recommend raising issues through the
`Github issues page <https://github.com/PMCC-BioinformaticsCore/janis/issues>`__ for Pipeline related issues.

Information about the project structure and more on contributing can be found within the documentation.


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

   tutorials/tutorial0
   tutorials/tutorial1
   tutorials/tutorial2
   tutorials/tutorial3
   tutorials/running
   tutorials/container

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Registry

   pipelines/index
   tools/index
   datatypes/index
   templates/index

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: API documentation

   references/index
   references/workflow
   references/commandtool
   references/tools/pythontool
   references/scattering
   references/secondaryfiles

.. toctree::
   :hidden:
   :maxdepth: 0
   :caption: Guides

   references/engines
   references/cwl
   references/wdl
   references/faq
   references/errors

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Development

   development/index
   development/contributing
   development/testing
   development/releasing

Indices and tables
===================


* :ref:`genindex`
* :ref:`search`
* :ref:`modindex`
