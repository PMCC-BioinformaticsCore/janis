.. WEHI Pipeline Definition documentation master file, created by
   sphinx-quickstart on Thu Jan 17 08:39:30 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Janis!
********************************

.. image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master
    :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis

.. image:: https://readthedocs.org/projects/janis/badge/?version=latest
    :target: https://janis.readthedocs.io/en/latest/?badge=latest)

.. image:: https://badge.fury.io/py/janis-pipelines.svg
    :target: https://pypi.org/project/janis-pipelines/

.. image:: https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/PMCC-BioinformaticsCore/janis

.. note::
	This repository is a work-in-progress. 	

Janis is an open source project to assist in writing and checking analysis pipelines at *transpile* time!

It was developed as part of the Portable Pipelines Project, a collaboration between:

- `Melbourne Bioinformatics (University of Melbourne) <https://www.melbournebioinformatics.org.au/>`_
- `Peter MacCallum Cancer Centre <https://www.petermac.org/>`_ 
- `Walter and Eliza Hall Institute of Medical Research (WEHI) <https://www.wehi.edu.au/>`_.

There are a few parts to this project that will help you to write high quality pipelines:

- **Pipeline** - A Python module that can be imported and provides a class structure to build a workflow.

- **Tool registry** - The tool registry contains tool definitions such as GATK4, BWA and other tool suites. These make up the *steps* in your workflow. Each tool has an exact docker container so you can be sure you're always using the correct tool. 

- **Types registry** - Sometimes a singular file isn't enough to contain all information; in some cases reference files have indexes or files as metadata, called `secondary files <https://docs.sevenbridges.com/page/secondary-files/>`_. To aid in sanity checking pipelines, on top of the ordinary data types (such as strings, files, numbers) we've created a number of collection types so you don't need to think about forgetting index files again!


Quick-start
===========

You can also Janis via PIP: ``pip install janis-pipelines``


Introduction
============

Before writing pipelines, it's useful to have some background knowledge of what makes up a Workflow. Simply put:

    *A workflow is a series of steps that are joined to each other.*

A diagram is great way to get an intuitive understanding of what a workflow is, for example to bake `Chocolate chip cookies <https://www.taste.com.au/recipes/chocolate-chip-cookies-2/1bfaa0e6-13b4-489d-bbd8-1cc5caf1fa32 />`_, we can refer to the following workflow:

.. image:: resources/bakecookie-wf.png

There are 8 steps to baking (and eating) a chocolate chip cookie, each represented by a box. The arrows represent the dependencies between the steps (what must be completed before the step can be started). For example, to beat the butter and the eggs, we must have acquired the ingredients. Or to place the rolled cookie balls onto the tray, we must have lined the 2 baking trays and mixed in the flour and chocolate.

We easily see how dependencies are depicted in a workflow, and given enough resources, some of these tasks could be completed at the same time. For example, we can preheat the oven while we're preparing the cookie dough.

We'll come back to this example a few times over the course of this introduction.


Inputs and Outputs
------------------

The inputs to a step is what we (the user) must provide to the worfklow in order to run it. This is 


Steps
-----


Reproducability
---------------


Portability
-----------

-------


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self

   workflowspecifications

   userguide

   buildingTools

   tools/index
   datatypes/index
   references/pipeline

Indices and tables
^^^^^^^^^^^^^^^^^^

* :ref:`genindex`
* :ref:`search`
.. * :ref:`modindex`
