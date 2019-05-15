# Running your workflow

Now you've built your workflow, let's run it on your computer.

> We're working on a runner for engine that allows you to run your favourite workflow engine without going through this manual export process. Stay tuned for more!

## Overview

In this tutorial, we're going to take the workflow that you've built using `janis`, export the workflow to a workflow specification language (CWL or WDL) and run it! 

#### Engines

To run the workflow, we're going to use a workflow execution engine, or just _engine_ for short. Some engines only accept certain workflow specifications, for example the CWL reference engine ([cwltool](https://github.com/common-workflow-language/cwltool)) only accepts CWL, while [Cromwell](https://github.com/broadinstitute/cromwell) (developed by the Broad Institute) accepts both WDL and CWL.

As you get the choice of which workflow specification you can export to, there are instructions for both on this guide.

### What you'll need

To complete this tutorial, you'll need to have an installation of janis (see the [Getting Started](/tutorials/gettingstarted) guide for more info). You can install janis by running:
```bash
pip3 install janis-pipelines
```

You'll also need a `janis` workflow (you can use the [`examples/simple`](/tutorials/simple) workflow), a workflow execution engine and a container software installed.


### Workflow Execution Engines

#### CWL
If you're using CWL, `cwltool` is a quick engine that is installable through `pip`. We'd recommend visiting the [install guide](https://github.com/common-workflow-language/cwltool#install), however you can install `cwltool` by running:

```bash
pip3 install cwltool
```

#### WDL

For WDL, we'd recommend `Cromwell`, developed by the Broad Institute and developed alongside the [openwdl](#) community. You can find their docs [here](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/). Cromwell is packaged as a `jar` file, you must have a Java 8 runtime installed (ie, running `java -version` in the console gives you `1.8.0` or higher).

You can download the Cromwell jar from their [release page](https://github.com/broadinstitute/cromwell/releases).

### Containerisation

Most commonly, people use [`docker`](https://www.docker.com/get-started) which both `cwltool` and `Cromwell` work with out of the box. Docker doesn't always work in all environments though, some HPC friendly alternatives include`singularity` and `udocker`.

However not all tools accept other containerisation technology as friendly, refer to the subsections for more information.

#### CWLTool

CWLTool does allow some user-space docker replacements out of the box, including both `singularity` and `udocker`, you can find instruction on usage on their [README (#using-user-space-replacements-for-docker)](https://github.com/common-workflow-language/cwltool#using-user-space-replacements-for-docker).

#### Cromwell

Cromwell is a little bit trickier, you'll need to build a backend configuration and override the `submit-docker` script. Cromwell have a guide on [getting started with containers](https://cromwell.readthedocs.io/en/stable/tutorials/Containers/) with examples for some batch systems. We'd really recommend you read through that page if you intend on using containers with Cromwell.


## Let's get started

We'll assume you have a workflow in the variable `w` (annotated with the `janis.Workflow` type for clarity):

```python
import janis

...

w: janis.Workflow = MyWorkflow()
```

### Exporting the translation

Let's write a CWL representation of our workflow to the console, we can do this by running:
```python
w.translate("cwl")
```

Awesome, you should see the workflow, then the tools and finally a blank inputs file written to your console! We can take a look at the method signature for `translate` to find out what else we can do:

```python
janis.workflow.translate(
    translation: SupportedTranslation,		# 'cwl' | 'wdl'
    to_console=True,						# Logs workflow to console
    to_disk=False,							# Write workflow + tools to disk
    export_path=ExportPathKeywords.default,	# Default: ~/Desktop/{name}/{language}/
    write_inputs_file=False,				# Write the inputs file to disk
    should_validate=False,					# Use $cwltool | $womtool to validate wf
    should_zip=True							# zip the tools/ into archive
    with_docker=True,						# Include docker containers export
    with_resource_overrides=False,			# generate resource overrides at a workflow level
)
```

#### Exporting to disk

To export a workflow in CWL for running, we can use the following statement:

> If you're using Cromwell or another engine that requires the tool dependencies as a zip, include the ` should_zip=True` argument.
```python
w.translate("cwl", to_disk=True, write_inputs_file=True)
```

After this runs, we'll have a few files (and a directory) within the export path:

- `${name}.cwl` - The workflow.
- `${name}-job.yml` - The inputs file.
- `tools/` - Directory of the dependencies.
- `tools.zip` - if you choose to include the `should_zip=True` argument.

The structure is identical for WDL, except with a `.wdl` extension.

### Using a workflow execution engine

Now that we have exported tool definitions, we can run our workflows using a workflow engine.

#### CWL

For CWL, change directory (`cd`) into where your workflow (and especially the `tools`) were exported, and then run the following command, replacing `${name}` with your workflow's name:

```bash
cwltool ${name}.cwl ${name}-job.yml
```

#### WDL
Similarly, you can run the following command to run a workflow in Cromwell.

```bash
java -jar $cromwell run ${name}.wdl -i ${name}-job.wdl -p tools.zip
```
You'll see the running of the workflow in the console, and the outputs 

Visit Cromwell's documentation for more ways to run workflows.

## Finished!

Congratulations on running your workflow! You can now get more advanced with the following tutorials:

- [Building a simple bioinformatics workflow](/tutorials/alignsortedbam)
- [Building tools](/tutorials/buildtools)

----

## Advanced arguments

### CPU and Memory overrides

- `runtime_cpu: int` is the number of cpus that are required.
- `runtime_memory: float` is specified in number of `Gigabytes`.

When exporting the workflow, you have the ability to generate schema that allows you to override the `cpu` and `memory` values at a workflow inputs level. That is, you can include parameters in your input file that will ask the engine to run your workflow with specific memory and cpu requirements.

You can generate a workflow with resource overrides like so:
```python
w.translate("cwl", with_resource_overrides=True, **other_kwargs)
```

You can generate the override inputs by executing the following method on the workflow:
```python
w.generate_resources_file(
    translation: translations.SupportedTranslation, # 'cwl' | 'wdl'
    hints: Dict[str, Any]=None,		# dictionary of hints
    to_console=True					# write the resources file to disk
)
```

In CWL, this looks like:
```yaml
taskname_runtime_cpu: null
taskname_runtime_memory: null
subworkflowname_task2name_runtime_cpu: null
subworkflowname_task2name_runtime_memory: null
```

And in WDL:
```json
{
    "$name.taskname_runtime_memory": null, 			
    "$name.taskname_runtime_cpu": null,
    "$name.taskname_runtime_disks": null,
    "$name.subworkflowname_task2name_runtime_cpu": null,
    "$name.subworkflowname_task2name_runtime_cpu": null,
    "$name.subworkflowname_task2name_runtime_disks": null
}
```

### Overriding export location
> Default: `export_path=~/Desktop/{name}/{language}/`

You can override the export path to a more convenient location, by providing the ``parameter to `translate`. Prefixing the path with `~` (_tilde_) will replace this per [`os.path.expanduser`](https://docs.python.org/3/library/os.path.html#os.path.expanduser). You can also use the following placeholders (see [`github/janis/translation/exportpath`](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/janis/translations/exportpath.py) for more information):

- `"{language}"` - 'cwl' or 'wdl'
- `"{name}"` - Workflow identifier
