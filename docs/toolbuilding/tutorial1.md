
# Tutorial 2 - Wrapping a new tool

## Introduction

A CommandTool is the interface between Janis and a program to be executed. Simply put, a CommandTool has a name, a command, inputs, outputs and a container to run in. Inputs and arguments can have a prefix and / or position, and this is used to construct the command line.

The Janis documentation for the [CommandTool](https://janis.readthedocs.io/en/latest/references/commandtool.html) gives an introduction to the tool structure and a template for constructing your own tool. A tool wrapper must provide all of the information to configure and run the specified tool, this includes the `base_command`, [janis.ToolInput](https://janis.readthedocs.io/en/latest/references/commandtool.html#tool-input), [janis.ToolOutput](https://janis.readthedocs.io/en/latest/references/commandtool.html#tool-output), a `container` and its version.

## Requirements

You must have Python 3.6 and Janis installed:

```bash
pip3 install janis-pipelines
```

You can check you have the correct version of Janis installed by running:

```bash
$ janis -v
--------------------  ------
janis-core            v0.7.1
janis-assistant       v0.7.8
janis-unix            v0.7.0
janis-bioinformatics  v0.7.1
--------------------  ------
```


### Container

> _Further information_: [Containerising your tools](https://janis.readthedocs.io/en/latest/tutorials/container.html)

> Guide on using containers

For portability, we require that you specify an OCI compliant `container` (eg: Docker) for your tool. Often they're will already be a container with some searching, however here's a guide on [preparing your tools in containers](https://janis.readthedocs.io/en/latest/tutorials/container.html) to ensure it works across all environments. 


## Samtools flagstat

In this workshop we're going to wrap the `samtools flagstat` tool. 

### Samtools project links

- Latest version: `1.9`
- Project page: [http://www.htslib.org/doc/samtools.html](http://www.htslib.org/doc/samtools.html)
- Github: [samtools/samtools](https://github.com/samtools/samtools)
- Docker containers: [quay.io/biocontainers/samtools](https://quay.io/repository/biocontainers/samtools?tag=latest&tab=tags) (automatically / community generated)
    - Latest tag: `1.9--h8571acd_11`


### Command to build

We want to replicate the following command for `Samtools Flagstat` in Janis:

```bash
samtools flagstat [--threads n] <in.bam>
```

Hence, we can isolate the following information:

- Base commands: `"samtools"`, `"flagstat"`
- The positional `<in.bam>` input
- The configuration `--threads` input


### Command tool template

The following template is the minimum amount of information required to wrap a tool. For more information, see the [CommandTool documentation](https://janis.readthedocs.io/en/latest/references/commandtool.html).

```python
from typing import List, Optional, Union
import janis as j

class ToolName(j.CommandTool):
    @staticmethod
    def tool() -> str:
        return "toolname"
        
    @staticmethod
    def base_command() -> Optional[Union[str, List[str]]]:
        pass
        
    @staticmethod
    def container() -> str:
        return ""

    @staticmethod
    def version() -> str:
        pass

    def inputs(self) -> List[j.ToolInput]:
        return []

    def outputs(self) -> List[j.ToolOutput]:
        return []
```

### Tool information

Let's start by creating a file with this template:

```bash
vim samtoolsflagstat.py
```

We can start by filling in the basic information:

- Rename the class to be `SamtoolsFlagstat`
- Add a unqiue tool identifier in the `tool()` method
- Fill in the `base_command` to be `["samtools", "flagstat"]`
- Add the container `"quay.io/biocontainers/samtools:1.9--h8571acd_11"`
- Add the version: `"1.9.0"`

You'll have a class definition like the following
```python
class SamtoolsFlagstat(j.CommandTool):
    @staticmethod
    def tool() -> str:
        return "samtoolsflagstat"

    @staticmethod
    def base_command() -> Optional[Union[str, List[str]]]:
        return ["samtools", "flagstat"]
    
    @staticmethod
    def container() -> str:
        return "quay.io/biocontainers/samtools:1.9--h8571acd_11"

    @staticmethod
    def version() -> str:
        return "1.9.0"
    
    # inputs and outputs
```

### Inputs


We'll use the [ToolInput](https://janis.readthedocs.io/en/latest/references/commandtool.html#tool-input) class to represent these inputs. A `ToolInput` provides a mechanism for binding this input onto the command line (eg: prefix, position, transformations). See the documentation for more ways to configure a ToolInput.

Our positional input is a Bam, so we'll import the Bam type from `janis` with the following line:

```python
from janis.data_types import Bam
```

Then we can declare our two inputs:

1. Positional bam input
2. Threads configuration input with the prefix `--threads`


```python
    # in the class
    def inputs(self) -> List[j.ToolInput]:
        return [
            # 1. Positional bam input
            j.ToolInput(
                "bam", 
                Bam, 
                position=1, 
                doc="Input bam to generate statistics for"
            ),
            # 2. `threads` inputs
            j.ToolInput(
                "threads", 
                j.Int(optional=True), 
                prefix="--threads", 
                doc="(-@)  Number of additional threads to use [0] "
            )
        ]
```


### Outputs

We'll use the [ToolOutput](https://janis.readthedocs.io/en/latest/references/commandtool.html#tool-output) class to collect and represent these outputs. A `ToolOutput` has a type, and if not using `stdout` we can provide a `glob` parameter.

The only output of `samtools flagstat` is the statistics that are written to `stdout`. We give this the name `"stats"`, and collect this with the `j.Stdout` data type:

```python
    # in the class
    def outputs(self) -> List[j.ToolOutput]:
        return [
            ToolOutput("stats", j.Stdout)
        ]
```


### Tool definition

Putting this all together, you should have the following tool definition:

```python
from typing import List, Optional, Union
import janis as j
from janis.data_types import Bam

class SamtoolsFlagstat(j.CommandTool):
    @staticmethod
    def tool() -> str:
        return "samtoolsflagstat"

    @staticmethod
    def base_command() -> Optional[Union[str, List[str]]]:
        return ["samtools", "flagstat"]

    @staticmethod
    def container() -> str:
        return "quay.io/biocontainers/samtools:1.9--h8571acd_11"

    @staticmethod
    def version() -> str:
        return "1.9.0"

    def inputs(self) -> List[j.ToolInput]:
        return [
        j.ToolInput(
             "bam", 
             Bam, 
             position=1, 
             doc="Input bam to generate statistics for"
        ),
        # 2. `threads` inputs
        j.ToolInput(
            "threads", 
            j.Int(optional=True), 
            prefix="--threads", 
            doc="(-@)  Number of additional threads to use [0] "
        )
    ]

    def outputs(self) -> List[j.ToolOutput]:
        return [j.ToolOutput("out", j.Stdout)]
```

## Testing the tool

We can test the translation of this from the CLI:

> If you have multiple command tools or workflows declared in the same file, you will need to provide the `--name` parameter with the name of your workflow.

```bash
janis translate samtoolsflagstat.py wdl # or cwl
```

In the following translation, we can see the WDL representation of our tool. In particular, the `command` block gives us an indication of how the command line might look:
```
task samtoolsflagstat {
  input {
    Int? runtime_cpu
    Int? runtime_memory
    File bam
    Int? threads
  }
  command {
    samtools flagstat \
      ${"--threads " + threads} \
      ${bam}
  }
  runtime {
    docker: "quay.io/biocontainers/samtools:1.9--h8571acd_11"
    cpu: if defined(runtime_cpu) then runtime_cpu else 1
    memory: if defined(runtime_memory) then "${runtime_memory}G" else "4G"
    preemptible: 2
  }
  output {
    File out = stdout()
  }
}
```



### Running the workflow

We can call the `janis run` functionality (default CWLTool) with our bam input:

```bash
janis run samtoolsflagstat.py --bam /data/brca1.bam
```

OUTPUT:
```
WID:        f9e89f
EngId:      f9e89f
Name:       samtoolsflagstatWf
Engine:     cwltool

Task Dir:   $HOME/janis/execution/samtoolsflagstatWf/20191114_155159_f9e89f/
Exec Dir:   None

Status:     Completed
Duration:   4s
Start:      2019-11-14T04:51:59.744526+00:00
Finish:     2019-11-14T04:52:03.869735+00:00
Updated:    Just now (2019-11-14T04:52:05+00:00)

Jobs: 
    [✓] samtoolsflagstat (N/A)       

Outputs:
    - out: $HOME/janis/execution/samtoolsflagstatWf/20191114_155159_f9e89f/output/out
2019-11-14T15:52:05 [INFO]: Exiting

```

We can test this by looking at our output file: 

```bash
cat $HOME/janis/execution/samtoolsflagstatWf/20191114_155159_f9e89f/output/out
```

```
20061 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
337 + 0 supplementary
0 + 0 duplicates
19971 + 0 mapped (99.55% : N/A)
19724 + 0 paired in sequencing
9862 + 0 read1
9862 + 0 read2
18606 + 0 properly paired (94.33% : N/A)
19544 + 0 with itself and mate mapped
90 + 0 singletons (0.46% : N/A)
860 + 0 with mate mapped to a different chr
691 + 0 with mate mapped to a different chr (mapQ>=5)
```

## Summary

- Learn about the structure of a CommandTool,
- Use an existing docker container,
- Wrapped the inputs, outputs and tool information in a Janis CommandTool wrapper,

### Next steps

- [Containerising a tool](https://janis.readthedocs.io/en/latest/tutorials/container.html) 
