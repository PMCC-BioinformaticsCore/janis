# Setup

In this step we'll install Janis on a local computer, setup a workflow engine and run a test workflow. 

> Instructions for setting up Janis on a compute cluster are under construction. 

## Requirements

- Local environment
- Python 3.6+
- Docker

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
    pip install janis-pipelines
    ```

3. Test that janis was installed:

    ```bash
	janis -v
	# --------------------  -------
    # janis-core            v0.7.3
    # janis-assistant       v0.7.10
    # janis-unix            v0.7.0
    # janis-bioinformatics  v0.7.1
    # --------------------  -------
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

You can test run an example workflow with Janis and CWLTool with the following command:

```bash
janis run --engine cwltool hello
```

You'll be presented with the progress screen as your workflow completes. Some things to note:

- `WID` - the janis identifier of your workflow.
- `Task Dir` - Where your workflow, output files and logs are.
- 

```
WID:        df5daa
EngId:      df5daa
Name:       hello
Engine:     cwltool

Task Dir:   $HOME/janis/hello/20191115_105042_df5daa/
Exec Dir:   None

Status:     Completed
Duration:   6s
Start:      2019-11-14T23:50:42.940196+00:00
Finish:     2019-11-14T23:50:48.453133+00:00
Updated:    Just now (2019-11-14T23:50:53+00:00)

Jobs: 
    [✓] hello (2s)       

Outputs:
    - out: $HOME/janis/execution/hello/20191115_105042_df5daa/output/out
```

There is a single output `out` from the workflow, cat-ing this result we get:

```bash
cat $HOME/janis/execution/hello/20191115_105042_df5daa/output/out
# Hello, World
```

### Overriding an input

The workflow `hello` has one input `inp`. We can override this input by passing `--inp $value` onto the end of our run statement, eg:

```bash
janis run --engine cwltool hello --inp "Hello, yourname"
# out: Hello, yourname
```


## Summary

- Setup a virtualenv
- Installed Janis and CWLTool
- Ran a small workflow with custom inputs

### Next steps

- [Workflow construction tutorial](https://janis.readthedocs.io/en/latest/tutorial1/construction.html)

