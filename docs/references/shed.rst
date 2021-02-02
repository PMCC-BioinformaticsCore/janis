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


Development Decisions
-----------------------

The design of the Janis tool repositories is a little strange, but we'll discuss some of the design decisions we made to give context.

We know that the underling command line for tools is pretty constant between versions, and we didn't want to create a new tool definition for each version for a few reasons:

- Easy to fix errors if code is in one place
- Makes pull requests easy when there's minimal code changes

For groups of tools that might share common inputs (eg: ``GATK4``), we wanted to declare a `Gatk4Base <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/blob/master/janis_bioinformatics/tools/gatk4/gatk4toolbase.py>`_, where we can limit the memory of each GATK tool by inserting the parameter ``--java-options '-Xmx{memory}G'``. Doing this in the GATK base means it's implemented in one location and ALL of the GATK tools get this functionality for free.


For these reasons, janis-bioinformatics uses a class inheritance structure, but this sometimes makes for some awkward imports.

1. Often the tool provider will declare some base, with some handy methods. For example, our Samtools base command is ``["samtools", SamtoolsCommand]``, so we'll create an abstract method for ``samtools_command`` that each tool must implement, and we can take care of the base command for them. We'll also declare a version of the base GATK with the version and containers method filled, we'll see why in a second.

.. code-block:: python

   import abc

   class SamToolsBase(abc.ABC):
       @abc.abstractmethod
       def samtools_command(self):
           pass

       def base_command(self):
           return ["samtools", self.samtools_command()]

       # other handy methods

   class SamTools_1_9:
       def container(self):
           return "quay.io/biocontainers/samtools:1.9--h8571acd_11"

       def version(self):
           return "1.9.0"


2. Now we can create our base definition for our class. Although a CommandTool usually requires a container and version, we'll use the one we declared earlier in the next step. Let's create a rough definition for Samtools Flagstat

.. code-block:: python

   class SamToolsFlagstatBase(SamToolsToolBase, ABC):
       def samtools_command(self):
           return "flagstat"

       def inputs(self):
           return [
               ToolInput("bam", Bam(), position=10)
           ]

       def outputs(self):
           return [ToolOutput("out", Stdout(TextFile))]

3. We can combine our Samtools version and FlagstatBase to create a version of flagstat, eg:

.. code-block:: python

   class SamToolsFlagstat_1_9(SamTools_1_9, SamToolsFlagstatBase):
       pass


Suggestions
++++++++++++++

See this GitHub issue for suggestions to improve the structure: `GitHub: janis-bioinformatics/issues/92 <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/issues/92>`_