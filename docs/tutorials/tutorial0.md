# Tutorial 0 - Introduction to Janis

Janis is workflow framework that uses Python to construct a declarative workflow. It has a simple workflow API within Python that you use to build your workflow. Janis converts your pipeline to the Common Workflow Language (CWL) and Workflow Description Language (WDL) for execution, and it’s also great for publishing and archiving.

Janis was designed with a few points in mind:

- Workflows should be easy to build,
- Workflows and tools must be easily shared (portable),
- Workflows should be able to execute on HPCs and cloud environments.
- Workflows should be reproducible and re-runnable.

Janis uses an *abstracted execution environment*, which removes the shared file system in favour of you specifiying all the files you need up front and passing them around as a File object. This allows the same workflow to be executable on your local machine, HPCs and cloud, and we let the `execution engine` handle moving our files. This also means that we can use file systems like ``S3``, ``GCS``, ``FTP`` and more without any changes to our workflow.

> Instructions for setting up Janis on a compute cluster are under construction. 

## Requirements

- Local environment
- Python 3.6+
- Docker
- Python virtualenv (`pip3 install virtualenv`)

> **NB**: This tutorial requires Docker to be installed and available.

## Installing Janis

We'll install Janis in a virtual environment as it preserves versioning of Janis in a reproducible way.

1. Create and activate the virtualenv:

    ```bash
    # create virtual env
	virtualenv -p python3 ~/janis/env
	# source the virtual env
	source ~/janis/env/bin/activate
    ```

2. Install Janis through PIP:

    ```bash
    pip install -q janis-pipelines
    ```

3. Test that janis was installed:

    ```bash
	janis -v
    # --------------------  ------
    # janis-core            v0.9.7
    # janis-assistant       v0.9.9
    # janis-unix            v0.9.0
    # janis-bioinformatics  v0.9.5
    # janis-pipelines       v0.9.2
    # janis-templates       v0.9.4
    # --------------------  ------
	```
	
### Installing CWLTool

[CWLTool](https://github.com/common-workflow-language/cwltool) is a reference workflow engine for the Common Workflow Language. Janis can run your workflow using CWLTool and collect the results. For more information about which engines Janis supports, visit the [Engine Support](https://janis.readthedocs.io/en/latest/references/engines.html) page.

```bash
pip install cwltool
```

Test that CWLTool has installed correctly with:

```bash
cwltool --version
# .../bin/cwltool 1.0.20190906054215
```


## Running an example workflow with Janis

First off, let's create a directory to store our janis workflows. This could be anywhere you want, but for now we'll put it at `$HOME/janis/`

```bash
mkdir ~/janis
cd ~/janis
```

You can test run an example workflow with Janis and CWLTool with the following command:

```bash
janis run --engine cwltool -o tutorial0 hello
```

You'll see the `INFO` statements from CWLTool in terminal.

> To see all logs, add `-d` to become: 
> ```bash
> janis -d run --engine cwltool -o tutorial0 hello
> ```

At the start, we see the two lines in our output:

```
2020-03-16T18:49:08 [INFO]: Starting task with id = 'd909df'
d909df
```

This is our workflow ID (wid) and is one way we can refer to our workflow.

After the workflow has completed (or in a different window), you can see the progress of this workflow with:

```bash
janis watch d909df

# WID:        d909df
# EngId:      d909df
# Name:       hello
# Engine:     cwltool
# 
# Task Dir:   $HOME/janis/tutorial0
# Exec Dir:   None
# 
# Status:     Completed
# Duration:   4s
# Start:      2020-03-16T07:49:08.367981+00:00
# Finish:     2020-03-16T07:49:11.881006+00:00
# Updated:    3h:51m:54s ago (2020-03-16T07:49:11+00:00)
# 
# Jobs: 
#     [✓] hello (1s)       
# 
# Outputs:
#     - out: $HOME/janis/tutorial0/out
```

There is a single output `out` from the workflow, cat-ing this result we get:

```bash
cat $HOME/janis/tutorial0/out
# Hello, World
```

### Overriding an input

The workflow `hello` has one input `inp`. We can override this input by passing `--inp $value` onto the end of our run statement. Note the structure for workflow parameters and parameter overriding:

```
janis run <run options> worklowname <workflow inputs>
```

We can run the following command:
```bash
janis run --engine cwltool -o tutorial0-override hello --inp "Hello, $(whoami)"

# out: Hello, mfranklin
```

### Running Janis in the background

You may want to run Janis in the background as it's own process. You could do this with `nohup [command] &`, however we can also run Janis with the `--background` flag and capture the workflow ID to watch, eg:

```bash
wid=$(janis run \
    --background --engine cwltool -o tutorial0-background \
    hello \
    --inp "Run in background")
janis watch $wid
```


## Summary

- Setup a virtualenv
- Installed Janis and CWLTool
- Ran a small workflow with custom inputs

### Next steps

- [Workflow construction tutorial](https://janis.readthedocs.io/en/latest/tutorials/tutorial1.html)

