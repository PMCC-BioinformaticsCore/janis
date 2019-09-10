Janis Tutorial - Introduction
#############################

Welcome to the tutorial for Janis! This tutorial takes the workshop material and presents it in a step-by-step process to learn Janis, and get your workflows running.

Janis is workflow framework that uses Python to construct a declarative workflow. It exposes a workflow API you use within Python to declare your workflow. Janis supports converting your pipeline to the Common Workflow Language (CWL) and Workflow Description Language (WDL) for execution and archival.

Janis was designed with a few priorities:

- Workflows and tools must be easily shared (portable).
- Execution must be able to occur on HPCs and cloud environments.
- Workflows should be reproducible and re-runnable.

Janis uses an *abstracted execution environment* which allows the workflow to be executable on your local machine, HPCs and cloud. This means, that instead of dealing with file paths, we just consider a file object and let the ``execution engine`` handle moving our files. This also means that we can use file systems like ``S3``, ``GCS``, ``FTP`` and more without any code changes.

.. note::
   At the end of the execution, we can use Janis to `move our output files <#>`_ into a predictable output structure.

This tutorial will introduce you to workflow concepts, and you'll construct a simple alignment workflow with Janis.

.. toctree::
   :maxdepth: 2

   setup
   construction

