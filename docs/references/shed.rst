Janis Shed and Toolboxes
==========================

One important component of Janis is the Toolbox, a set of command tools and workflows that are immediately available to you. There are toolboxes for each domain, and we refer to the collection of these toolboxes as a "JanisShed".

In essence, they're Python packages that Janis knows about and can use, currently there are two toolboxes:

- Unix (`GitHub <https://github.com/PMCC-BioinformaticsCore/janis-unix>`__ | `Docs <https://janis.readthedocs.io/en/latest/tools/unix/index.html>`__)
- Bioinformatics (`GitHub <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`__ | `Docs <https://janis.readthedocs.io/en/latest/tools/bioinformatics/index.html>`__)


Adding a tool
--------------

To add a tool to an existing toolbox, please refer to `Adding a tool to a toolbox <https://janis.readthedocs.io/en/latest/tutorials/toolbox.html>`_.


Custom Toolbox
---------------

You might have a set of tools that you you want to use in Janis (and take advantage of the other features), but may not be allowed to publish these (and they're container). In this case, we'd recommend you set up your own toolbox as a Python module, and by setting a few entrypoints Janis will be able to pick up your tools.

We won't discuss how to setup a Python package from scratch, but you might want to follow the existing toolbox repos as a guide. Janis uses entrypoints to find existing toolboxes, here's a great guide on `Python Entry Points Explained <https://amir.rachum.com/blog/2017/07/28/python-entry-points/>`__.

Specifically, you'll need to add the lines to your ``setup.py`` (replacing ``bioinformatics`` with a keyword to describe your extension, and ``janis_bioinformatics`` with your package name):

.. code-block:: none

   entry_points={
       "janis.extension": ["bioinformatics=janis_bioinformatics"],
       "janis.tools": ["bioinformatics=janis_bioinformatics.tools"],
       "janis.types": ["bioinformatics=janis_bioinformatics.data_types"],
   },

Currently, the Janis assistant doesn't support configuring Cromwell or CWLTool to interact with private container repositories. If this is something you need, I'd recommend reaching out on GitHub and Gitter, or translating your output to run the workflow on your own engine.

