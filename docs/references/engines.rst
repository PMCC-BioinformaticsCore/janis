Engine support
==========================

Janis currently support 2 engines:

- CWLTool
- Cromwell


Both of these engines have their limitations, and there will instructions in this document about those limitations and how they can be installed / configured.


CWLTool
--------

`CWLTool <https://github.com/common-workflow-language/cwltool>`_ is the reference implementation for the Common Workflow Language (CWL). It's command line driven and can use Docker / Singularity for container driven processes. Unsurprisingly Janis transpiles your workflow to CWL to run with CWLTool.


Limitations
++++++++++++

CWLTool only runs one tool at a time in serial. Even though there is a *parallel* mode, this isn't production stable, and Janis can't provide accurate metadata and progress tracking with this mode. CWLTool by itself also doesn't have the ability to submit to batch systems. With some configuration, the manager of Janis (that spins up CWLTool) can be submitted to a Batch System.

Janis determines the progress of your workflow by processing the stdout and converting it into the internal Workflow metadata models.

Known issues:

- `Scattering on multiple fields with secondary files <https://github.com/common-workflow-language/cwltool/issues/1208>`_


Installation notes
++++++++++++++++++

The base CWLTool can be installed through Pip, see other installation methods here: `GitHub: common-workflow-language.cwltool <https://github.com/common-workflow-language/cwltool#install>`_.

In addition **Node** should be installed, otherwise CWLTool is going to try and call out to Docker to run a node server. This doesn't always work in shared environments.


Cromwell
---------

Cromwell is a workflow engine produced by the Broad Institute. It is very customisable for shared file systems (HPCs) and supports cloud computation through GCP (and supposed AWS).

Limitations:

CWL provides no pub-sub method for workflow notifications. Progress tracking is provided through requesting the metadata from the REST endpoint, and transforming this metadata into a recognised format for Janis. This metadata becomes quite large as the workflow grows, so is polled every 5-6 seconds.

- `CWL: InitialWorkDirRequirement not returning new localized path when constructing command <https://github.com/broadinstitute/cromwell/issues/4775>`_
- `CWL output doesn’t multiple workflow outputs that reference the same output of a tool <https://broadworkbench.atlassian.net/browse/BA-5790>`_


Installation notes:
+++++++++++++++++++

Cromwell requires a Java 8 runtime environment, this must be loaded at submit for Janis to start Cromwell. Janis looks in the following places for a Cromwell jar:

- Environment (path): ``JANIS_CROMWELLJAR``
- Configuration (path) : ``cromwell.jarpath``
- ConfigDir (globbed): Searches inside the ``configDir`` for the pattern: ``"cromwell-*.jar"`` (default ``~/.janis/``).

If no version of cromwell is found, then the latest version is downloaded from GitHub and placed in the config dir.

MySQL
*******

Additionally, Janis (*expected* >= v0.8.0) will attempt to start a MySQL container for persistence, stability and eventually call caching of Cromwell. This is facilitated using Docker or Singularity. By default this is Docker, and can currently only be configured using an environment template.

The data files for MySQL are stored within your task execution directory, as a workflow gets larger these database files might grow in size. MySQL can be disabled using the ``--no-mysql`` flag on the CLI for ``run``.


Unsupported engines
-------------------

Here are a couple of other engines that we tried, but either haven't got them working yet or they don't quite fit our purposes.

Toil
++++++

Toil is primarily Python API driven workflow, it supports CWL through a bridge that uses CWLTool to parse and run the jobs. Although it is a good candidate for Janis, it's difficult to determine the progress of a workflow through the internal job store. It might be possible to generate the workflow in a particular way and rewrite the CWLToil bridge and parse the --stats output, but it was infeasible for this project.

- `How to figure out progress using --stats <https://github.com/DataBiosphere/toil/issues/2580>`_
- `toil-cwl-runner isn't assigning defaults for subworkflow <https://github.com/DataBiosphere/toil/issues/2727>`_.


MiniWDL
++++++++

MiniWDL is an alternative engine for WDL support. miniWDL runs the exemplar WGS pipelines as expected on targeted panels, but no investigation has been completed how to obtain progress information or whether it can be configured to work with batch systems, the cloud or with singularity.




