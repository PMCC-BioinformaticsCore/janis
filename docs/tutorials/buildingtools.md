  
# Building Tools   

_This page is under construction_  

Welcome to this tutorial on building a a tool. Prior to starting this, we recommend completing:

- [Janis Tutorial](/tutorial1/)

We also recommend reading the following documentation page on Command Tools: 

- [Command Tool API reference](../references/commandtool.html)

## Overview

Janis was designed to be backed by a strongly typed tool registry that is made easily available. This set of documentation includes [automatically generated tool pages](/tools/). 

As you're writing workflows, you'll find a tool that you want to place in your workflow, but it won't be in our toolshed. In this tutorial we're going to analyse a tool to determine its inputs and outputs to build a tool wrapper that everyone can use.

The tools in Janis take a bit of effort to build, but they're designed to be built once and reused as much as possible.

### What you'll need

You'll need a workfing knowledge of the tool and a container. 

#### Container

> _Further information_: [Dockerizing your tools](https://janis.readthedocs.io/en/latest//tutorials/docker.html)

For portability, we require that you specify a `Docker` (or other OCI compliant) container for your tool. Often they're will already be a container with some searching, however here's a guide on [preparing your tools as dockers](https://janis.readthedocs.io/en/latest/tutorials/docker.html) to ensure it works across all environments and workflow specifications.

Your container must be available to your execution engine and as such we'd recommend hosting your container on public repository such as Dockerhub / quay.io. 

## Let's get started!

Janis is built with extensibility in mind. We currently have two registries of prebuilt tools: 

- [`janis-bioinformatics`](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics)
- [`janis-unix`](https://github.com/PMCC-BioinformaticsCore/janis-unix)

It's great to contribute your tool into these repositories to avoid repetitive rebuilding of tools, however this isn't required. Your tool only has to be available to your workflow (through imports or in the same file).
 
### Preparation
  
The recommended structure of a tool is the following:  
```  
tool_repo/  
├── __init__.py
├── tool_producer/  
│   ├── __init__.py
│   ├── tool1_name/  
│   │   ├── base.py  
│   │   ├── versions.py
│   └── tool2_name/  
│   │   ├── base.py  
│   │   ├── versions.py
│   ... other tools  
```  

> The `__init__.py` files are important to ensure that your tool files are correctly built and distributed, we'll also use them later to export our new tool.

### Setting up the file:

A new tool definition must subclass the janis.CommandTool class and implement the required abstract methods:

> Based on [this template](../references/commandtool.html#template)

```python
from abc import ABC
from typing import List, Optional, Union
import janis as j

class ToolNameBase(j.CommandTool, ABC):
    @staticmethod
    def tool() -> str:
        return "toolname"

    @staticmethod
    def base_command() -> Optional[Union[str, List[str]]]:
        pass

    def inputs(self) -> List[j.ToolInput]:
        return []

    def outputs(self) -> List[j.ToolOutput]:
        return []


class ToolName_Version(ToolNameBase):

    @staticmethod
    def container() -> str:
        return ""

    @staticmethod
    def version() -> str:
        pass
```

In addition, we can include the following metadata:

```python
    # within class 
    def arguments(self) -> List[j.ToolArgument]:
        # parameters that are not overridable
        return []

    def friendly_name(self) -> str:
        pass

    @staticmethod
    def tool_module() -> str:
        pass

    @staticmethod
    def tool_provider() -> str:
        pass
        
    def bind_metadata(self) -> j.ToolMetadata:
       pass
```

### Tool inputs

A tool input is a named input to a tool that has four primary attributes:
- `tag: str` – The identifier of the input (unique to inputs and outputs)
- `input_type: ParseableType` - The data type that this input accepts
- `position: int` - The position of the input to be applied. (Default = 0)
- `prefix: string` - The prefix to be appended before the element. (By default, a space will also be applied, see separate_value_from_prefix for more information)

_You must provide a position or prefix for an input to be bound onto the command line._

There are ways to customise how a `ToolInput` is bound onto the command line:
- `separate_value_from_prefix` (default: `True`) - eg: set to false when you expect the command line to look like `--prefix=myvalue`. 
- `prefix_applies_to_all_elements` – Applies the prefix to each element of the array (Array inputs only)
- `shell_quote` – Stops shell quotes from being applied in all circumstances, useful when joining multiple commands together.
- `separator` – The separator between each element of an array (defaults to `' '`)
- `localise_file` – Ensures that the file(s) are localised into the execution directory.
- `default` – The default value to be applied if the input is not defined.
- `doc` - Documentation string used for  

Example:
```python
ToolInput(
    tag="tumorBams",
    input_type=Array(BamBai),
    prefix="-I",
    prefix_applies_to_all_elements=True,
)
```

### Tool outputs

A tool output is a named output of a tool. The ToolOutput object instructs how the engine should collect an output and how it may be referenced in a workflow.

- `tag: str` – The identifier of a output, must be unique in the inputs and outputs.
- `output_type: ParseableType` – The type of output that is being collected.
- `glob: Selector` – How to collect this output, can accept any janis.Selector.

Example:

```python
ToolOutput(
    "stats",
    TextFile(extension=".stats"),
    glob=InputSelector("outputFilename") + ".stats",
    doc="To determine type",
),


```

## Resources (CPU / Memory)

> Further information: [_Resource Overrides_](/resourceoverrides.md), and the `hints` section of this page

Resources are an important part of optimising tools, and `janis` exposes a few ways to allow your tool to adapt with different inputs.

> Expressions are NOT supported yet.

There are 2 ways you might need to reference resources in your tool definition:

1. Your tool requires the CPU or Memory runtime values as parameters
2. To generate efficient resource files.

### Runtime values as Parameters

> Further information: [Runtime Parameters](

If your tool requires the CPU or Memory runtime values (maybe as Java args or as a thread parameter), you can provide the `CpuSelector` or `MemorySelector` classes as defaults. You instantiate them with no parameters (or a default), and they get the value provided by the `runtime_cpu` or `runtime_memory` parameters.

For example, in the [`bwa mem`](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/blob/master/janis_bioinformatics/tools/bwa/mem/base.py#L111) tool wrapper, the `thread` parameter has the following definition:

```python
ToolInput("threads", Int(optional=True), default=CpuSelector(), prefix="-t", doc="Number of threads. (default = 1)")
```

Which results in the following WDL:
```wdl
Int? threads = runtime_cpu
```
Or CWL:
```cwl
inputs:
- doc: Number of threads. (default = 1)
  id: threads
  inputBinding:
    prefix: -t
    valueFrom: $(inputs.runtime_cpu)
  label: threads
  type:
  - int
  - 'null'
```

### Generating a runtime override file
  
Resources are specified with the tools, but kept separate from the workflow specification. This allows us to generate a `runtime.yml` or `runtime.json` inputs file that  specifies all the resources you need to run the workflow.  When you build the tool, you might know exactly how the tool should react to different `hints`, and for this reason you can override the following method stubs:  
  
1. Memory: The memory must be specified in Gigabytes:  
  
```python  
def memory(self, hints: Dict[str, Any]) -> Optional[float]:  
    return None
 ```
    
2. CPUs: The number of CPUs must be specified in whole numbers  
  
```python  
def cpus(self, hints: Dict[str, Any]) -> Optional[int]:  
    return None
 ```  
  
## Advanced versioning  
_The section is under construction_.  

  
## Piped and combining commands  
  
When you want two tools to pipe data between them, unfortunately there's not a great solution. CWL allows you to tag inputs and outputs as [`streamable`](https://www.commonwl.org/v1.0/CommandLineTool.html#CommandInputParameter), however there are [no workflow engines](https://github.com/broadinstitute/cromwell/issues/3454#issuecomment-455367417) that properly support this.  
  
For this reason, the current (less than impressive) suggestion is to manually pipe these together by adding a `|` (pipe) `ToolArgument` at certain positions, and ensure that the positions of the other `ToolInput`'s are explicit and sequential.  
  
You will also need to add a `shell_quote=False` to the input of every `ToolInput` / `ToolArgument`, to ensure that when the command is generated by the CWL / WDL generators, there are no shell quotes added.  
  
For example:  

```python  
 Piped tool command: #    tool1_baseCommand -inp1 ${tl1_input} | tool2_baseCommand -inp2 ${tl2_input}  
  
class PipedTool(CommandTool):  
	def inputs(self):  
		 return [
			 ToolInput("tl1_input", String(), position=2, prefix="-inp1" shell_quote=False),  
			 ToolInput("tl2_input", String(), position=5, prefix="-inp2", shell_quote=False)  
		 ]

	def arugments(self):  
		return [
			ToolArgument("tool1_baseCommand", position=1, shell_quote=False),  
			ToolArgument("|", position=3, shell_quote=False),  
			ToolArgument("tool2_baseCommand", position=4, shell_quote=False)  
		 ] 
	# outputs and other required commands here
 ```  
### Using an input more than once

Sometimes you'll need to use an input more than one. An easy way to perform this is to declare a `ToolArgument` using a `InputSelector` as the value, selecting your input you want to redeclare. You are still able to declare a new position and prefix.
  
## Building a workflow as a tool  

Janis supports the embedding of subworkflows as tools. To ensure the documentation generator will automatically pick up the tool, we'll subclass the `Workflow` class, and place our inputs, steps, outputs and connections within the `__init__` method. And then all the other steps are the same in regards to exporting your workflow correctly.
 