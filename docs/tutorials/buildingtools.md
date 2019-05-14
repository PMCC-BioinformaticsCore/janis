  
# Building Tools   

_This page is under construction_  

Welcome to this tutorial on building a a tool. Prior to starting this, we recommend completing these tutorials:

- [Getting started with Janis](/tutorials/gettingstarted)
- [Constructing your first workflow](/tutorials/simple)

## Overview  

At some point, you'll find a tool that you want to place in your workflow, but it won't be in our toolshed. In this tutorial we're going to analyse a tool to determine its inputs and outputs to build a tool wrapper that everyone can use.

### What you'll need

To complete this tutorial, you'll need to have an installation of janis (see the [Getting Started](/tutorials/gettingstarted.html) guide for more info). You can install janis by running:
```bash
pip install janis-pipelines
```

#### Container

> _Further information_: [Dockerizing your tools](https://janis.readthedocs.io/en/latest//tutorials/docker.html)

For portability, we require that you specify a `Docker` (or other OCI compliant) container for your tool. Often they're will already be a container with some searching, however here's a guide on [preparing your tools as dockers](https://janis.readthedocs.io/en/latest//tutorials/docker.html) to ensure it works across all environments and workflow specifications.

## Let's get started!

### Preparation
  First off, you need to decide where the tool should sit. If it's a unix tool, it should sit under the [`janis`](https://github.com/PMCC-BioinformaticsCore/janis) repository, if it's a bioinformatics tool: under the [`janis-bioinformatics`](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics) repo. If it's a new domain, raise a [Github issue](https://github.com/PMCC-BioinformaticsCore/janis/issues/new).
  
The general structure of the tools are:  
```  
janis_bioinformatics/  
├── __init__.py
├── tool_producer/  
│   ├── __init__.py
│   ├── tool1_name/  
│   │   ├── toolname_version1.py  
│   └── tool2_name/  
│   │   ├── toolname_version2.py  
│   ... other tools  
```  

So, place your file called: `$toolname_$version.py`in a subfolder of the tool name, in a subfolder of the manufacturer / producer of the tool.

> The `__init__.py` files are important to ensure that your tool files are correctly built and distributed, we'll also use them later to export our new tool.

### Setting up the file:


 
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

See the [Building a simple bioinformatics workflow](/tutorials/alignsortedbam) tutorial as an example of building a workflow as a reference-able tool.
 
## Regenerating documentation