# Running your workflow

Now you've built your workflow, let's run it on your computer.

_We're working on a runner for engine that allows you to run your favourite workflow engine without going through this manual export process. Stay tuned for more!

## Overview

In this tutorial, we're going to run the workflow you have built with `janis-runner`. Behind the scenes this uses Cromwell or CWLTool to manage the execution.

#### Engines

To run the workflow, we're going to use a workflow execution engine, or just _engine_ for short. Some engines only accept certain workflow specifications, for example the CWL reference engine ([cwltool](https://github.com/common-workflow-language/cwltool)) only accepts CWL, while [Cromwell](https://github.com/broadinstitute/cromwell) (developed by the Broad Institute) accepts both WDL and CWL.

## What you'll need

To complete this tutorial, you'll need to have an installation of janis (see the [Getting Started](/tutorials/gettingstarted) guide for more info), a workflow to run and an engine available.

### Execution Engine

#### CWLTool

`cwltool` is the reference engine for CWL. We'd recommend visiting the [install guide](https://github.com/common-workflow-language/cwltool#install), however you can install `cwltool` through Pip by running:

```bash
pip3 install cwltool
```

#### WDL

For WDL, we'd recommend `Cromwell`, developed by the Broad Institute and developed alongside the [openwdl](#) community. You can find their docs [here](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/). Cromwell is packaged as a `jar` file, you must have a Java 8 runtime installed (ie, running `java -version` in the console gives you `1.8.0` or higher).

Janis will automatically download Cromwell, or it will automatically detect versions in `~/.janis/`. It's possible to configure Cromwell to run in a variety of ways through the [configuration guide](../runner/index.html).

### Containerisation

Most commonly, people use [`docker`](https://www.docker.com/get-started) which both `cwltool` and `Cromwell` support by default. Some HPC friendly alternatives include`singularity` and `udocker`; although these engines can be configured independently with both of these user-space replacements, we're still working on an easy way to configure Janis to use these.

## Let's get started

We're going to use the Janis command-line interface (CLI) to run these workflows and track their progress. For this example, we're going to run the [Hello](#) workflow and put the relevant output files in an output directory called `hello_dir`. 


By default, Janis will run the workflow using CWLTool, this can be changed by passing the `--engine` parameter with either `"cwltool"` or `"cromwell"`.

Let's do that now:

```bash
janis run --engine [cwltool|cromwell] -o hello_dir hello 
```

We'll see some important information in our console:

- The ID of our task. This is useful to gain access to our workflow information at a later date. 
```
2019-09-26T00:25:00+00:00 [INFO]: Starting task with id = 'a4f047'
```

- The metadata screen, this tells us about the progress of our workflow.
    - There are 3 status indicators:
        - `[...]`: processing
        - `[~]`: running
        - `[✓]`: completed
        - `[!]`: failed
```
TID:        a4f047
WID:        1ca1d418-df9b-45ea-b6b8-2d210be310b3
Name:       hello

Engine:     cromwell
Engine url: localhost:56232

Task Dir:   /Users/franklinmichael/janis/execution/hello/20190926_102526_a4f047/
Exec Dir:   /Users/franklinmichael/cromwell-executions/hello/1ca1d418-df9b-45ea-b6b8-2d210be310b3

Status:     Completed
Duration:   17
Start:      2019-09-26T00:25:55.771000+00:00
Finish:     2019-09-26T00:26:13.079000+00:00

Jobs: 
    [✓] hello.hello (13s :: 76.5 %)
       
...

Finished managing task 'a4f047'. View the task outputs: file:///Users/franklinmichael/janis/execution/hello/20190926_102526_a4f047/
a4f047
```

If we look at the contents of our output directory (`$task dir + "/outputs/"`, eg: `/Users/franklinmichael/janis/execution/hello/20190926_102526_a4f047/outputs/`), we see one file called `hello.out` which contains the string "Hello, World!".

### Overriding workflow inputs

The `hello` workflow allows us to override the string that gets printed to the console. We can see in the [`hello` documentation](https://janis.readthedocs.io/en/latest/tools/unix/hello.html), there is an input called `inp` (optional string) that we can override.

There are two ways to override inputs in Janis:

1. Creating a yaml job file: `janis inputs hello > hello-job.yml`, this file can be modified and provided to the run with `--inputs hello-job.yml`.
2. Parameter overrides within the run command.

We're going to take the second route by providing `"Hello, Janis!"` to the workflow:

```bash
janis run hello --inp "Hello, Janis!"
```

Outputting the result in the outputs folder yields:
```
Hello, Janis!
```

Success!





## Finished!

Congratulations on running your workflow! You can now get more advanced with the following tutorials:

- [Building a simple bioinformatics workflow](/tutorials/alignsortedbam)
- [Building tools](/tutorials/buildtools)

----

## Advanced arguments

### CPU and Memory overrides

- `runtime_cpu: int` is the number of cpus that are required.
- `runtime_memory: float` is specified in number of `Gigabytes`.
- `runtime_disks: string` is specified in the format `target size type` (eg: `"local-disk 100 SSD"`)

When exporting the workflow, you have the ability to generate schema that allows you to override the `cpu` and `memory` values at a workflow inputs level. That is, you can include parameters in your input file that will ask the engine to run your workflow with specific memory and cpu requirements.

You can generate a workflow with resource overrides like so:
```python
w.translate("cwl", with_resource_overrides=True, **other_kwargs)
```

You can generate an inputs file with the overrides withe following:
```python
def generate_inputs_override(
    additional_inputs=None, 
    with_resource_overrides=False, 
    hints=None
):
    pass
```

or with the CLI:

```bash
janis inputs myworkflow.py [--resources] > job.yaml
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
> Default: `export_path='.'`

You can override the export path to a more convenient location, by providing the ``parameter to `translate`. Prefixing the path with `~` (_tilde_) will replace this per [`os.path.expanduser`](https://docs.python.org/3/library/os.path.html#os.path.expanduser). You can also use the following placeholders (see [`github/janis/translation/exportpath`](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/janis/translations/exportpath.py) for more information):

- `"{language}"` - 'cwl' or 'wdl'
- `"{name}"` - Workflow identifier

You can output a workflow through the CLI with the following commmand:
```bash
janis translate myworkflow.py [cwl|wdl] --output-dir /path/to/outputdir/
```
