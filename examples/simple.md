
# Building a simple workflow

Let's start with something small, something simple!

## Overview 

> A completed version of this workflow is available within the janis repo: [`janis/examples/simple.py`](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/janis/examples/simple.py) 

 We have a _tar_ archive that contains a Java file. In this tutorial we're going to untar the archive and compile the result of the archive.

### You will need

A basic installation of janis is required (see the [Getting Started](https://janis.readthedocs.io/en/latest/tutorials/tutorial0.html) guide for more info). You can install janis by running:

```bash
pip3 install janis-pipelines
```

## Let's get started!

We've discussed what our workflow will do, let's identify the tools we'll use to build this. All of these unix tools are bundled with the core installation of `janis`.

- [`untar`](https://janis.readthedocs.io/en/latest/tools/unix/tar.html)
- [`cat`](https://janis.readthedocs.io/en/latest/tools/unix/cat.html)

These tools are linked to their generated documentation on this site. It's highly recommended you become familiar with this documentation website, as it is the easiest way to determine the tool's requirements.

### Setting up
 
Here's a template for our simple workflow.

Let's put this into a file called: `simple.py`.

```python
import janis

 # tool and data_type imports go here

class Simple(janis.Workflow):
  def id(self):
    return "simpleWorkflowIdentifier"

  def constructor(self):
    # workflow construction here
	
```

### Inputs and Outputs

We only have one input, a TarFile. We'll use the type TarFile (from `janis.data_types`). We can use the import statement:

```python
from janis.data_types import TarFile
```

And the declaration of our input:
```python
def constructor(self):
  # workflow construction here
  self.input("tarfile", TarFile)
```

### Steps

We can import the tools we're going to use from the unix toolshed:
```python
from janis.tools import Compile, Untar
```

We'll instantiate them and provide the inputs as the previous steps or the inputs. Some code here will make more sense:

```python
  def constructor(self):
    # workflow construction here
    self.input("tarfile", TarFile)
    
    # Untar has an input called 'tarfile', which we'll connect the input 'tarfile' to.
    self.step(
      "untar", 
      Untar(tarfile=self.tarfile)   # Use the input 'tarfile`
    )

    # - The input to Compile 'file' will be the result from 'untar.out'
    # - As the output of untar is an array and compile only takes one input
    #     we need to scatter the task for each 'file'.
    self.step(
      "compile", 
      Compile(file=self.untar.out), # use the result from 'untar'
      scatter="file"                # Scatter on 'file'
    )
```


### Outputs

We know we're going to want to output the result of the Tar, hence we only need one output. Janis will automatically determine the output type, you can simply add the following line beneath the previous step declarations:

```python
  def constructor(self):
    # workflow construction here

    # ... 
    self.output("compiled", source=self.compile)
```


### The final result

```python
import janis
# Data types - These help us logically connect workflows
from janis.data_types import TarFile
# Tools - The command line tools we're going to call
from janis.tools import Compile, Untar


class SimpleWorkflow(janis.Workflow):
    def id(self) -> str:
        return "simple"

    def friendly_name(self):
        return "Simple workflow: Untar, Compile, Tar"

    def constructor(self):
        self.input("tarfile", TarFile)
        self.step("untar", Untar(tarfile=self.tarfile))
        self.step("compile", Compile(file=self.untar.out), scatter="file")
        self.output("out", source=self.compile)
```

### Export

Now that we have a workflow, it's time to export it to a representation. The WDL translation gives
us a good idea of how the tools link together. We can translate this using the CLI for Janis:

```bash
janis translate simple.py wdl
```

Here's the exported workflow (no tools):

```wdl
version development

import "tools/untar.wdl" as U
import "tools/javacompiler.wdl" as J
import "tools/tar.wdl" as T

workflow simpleWorkflow {
  input {
    File tarFile
  }
  call U.untar as untar {
    input:
      tarFile=tarFile
  }
  # This is the scatter that we discussed between untar and tar
  scatter (f in untar.files) {
     call J.javacompiler as compile {
      input:
        file=f
    }
  }
  call T.tar as tar {
    input:
      files=[untar.files, compile.compiled]
  }
  output {
    File out = tar.tarred
  }
}
```

If you want to export the translation to disk including all the tools, you can do so with the `-o DIRECTORY` command, eg to output the translation in the current directory you could use:

```
janis translate simple.py wdl -o .
```
## Finished!

Congratulations on constructing and exporting your first workflow, and i You can now get more advanced with the following tutorials:

- [Building a simple bioinformatics workflow](https://janis.readthedocs.io/en/latest/tutorials/tutorial1.html)
- [Running workflows](https://janis.readthedocs.io/en/latest/tutorials/running.html)
