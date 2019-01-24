# python-wdlgen

[Workflow Description Language](http://www.openwdl.org) is way to describe tasks and workflows in a "_human readable and writable way_". It was initially developed and offered by [Broad Institute](https://software.broadinstitute.org/) to be paired with their workflow engine [Cromwell](https://cromwell.readthedocs.io/en/stable/), however it has since been made open source with other engines such as [Toil](http://toil.readthedocs.io/en/3.14.0/running/wdl.html) and [DNAnexus\*](https://github.com/dnanexus/dxWDL).

It's based almost exclusively off the version [1.0 Workflow Description Language specification](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#task-inputs).

## Motiviation

I needed an easy way to generate some _BASIC_ WDL through some in memory objects, and I was using ([a fork](https://github.com/illusional/python-cwlgen) of) [common-workflow-language/python-cwlgen](https://github.com/common-workflow-language/python-cwlgen), I figured I could open this up to see what use it has.

## Installation

**Goal:** Put on PIP or similar.

Otherwise gitsubmodules might be your saviour here unfortunately.

## General support

This software is provided as-is, without warranty of any kind ... and so on.

It's a pretty dumb wrapper that uses string interpolation to generate the structure. It wouldn't handle automatically escaping illegal characters.

Generally it supports:

- Types - All types are represented as a `WdlType`, which can either be a [`PrimitiveType`](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#types), or an `ArrayType` (see [goal](#goals)). Also supports the postfix quantifiers.

- Workflow creation (`wdlgen.Workflow`)
	- manual imports (`wdlgen.Workflow.WorkflowImport`)
	- inputs (`wdlgen.Input`)
	- outputs (`wdlgen.Output`)
	- calls:
		- general call (`wdlgen.WorkflowCall`)
		- scatter (`wdlgen.WorkflowScatter(WorkflowCall[])`)

- Task creation (`wdlgen.Task`) - This is based similar to how [CWL constructs its commands](https://www.commonwl.org/v1.0/CommandLineTool.html#CommandLineTool).
	- inputs: `wdlgen.Input`
	- outputs: `wdlgen.Output`
	- runtime: `wdlgen.Task.Runtime`
	- command: `wdlgen.Task.Command`
		- arguments: `wdlgen.Task.Command.Argument`
		- inputs: `wdlgen.Task.Command.Input`

## How to use

This will give you a _brief_ overview on how to _use_ python-wdlgen. Goals are to improve the write a proper documentation spec, but if you have a moderate understanding of workflows in either CWL or WDL, this code will hopefully be fairly intuitive.

Every class inherits from a `WDLBase` which means it must have a `get_string()` method which returns the string representation of the class, calling this on any children it may have.

### Types

All types are represented as a WDLType, which has a parse method. It's a little overkill in some cases, but makes managing attributes a bit easier.

```python
parsed_string = wdlgen.WdlType.parse("String")	# WdlType<PrimitiveType<String>>
parsed_op_str = wdlgen.WdlType.parse("String?") # WdlType<PrimtiveType<String>>
parsed_array = wdlgen.WDLType.parse("File[]")	# WdlType<ArrayType<File>>
parsed_ar_oq = wdlgen.WdlType(parse("Int?[]+"))	# WdlType<ArrayType<Int?> (+)>
```

You can also construct these manually:
```python
parsed_string = WdlType(PrimitiveType("String"))
parsed_op_str = WdlType(PrimtiveType("String", optional=True))
parsed_array = WdlType(ArrayType(WdlType(PrimitiveType("File"))))
parsed_ar_q = WdlType(ArrayType(WdlType(PrimitiveType("Int"), optional=True), requires_multiple=True))
```

### Input / Output
Input: `wdlgen.Input(data_type: WdlType, name: str, expression: str = None)`

Output: `wdlgen.Output(data_type: WdlType, name: str, expression: str = None)`

both of which output something like:
	`{WdlType} {name} [= {expression}]`

### Task

A task is a collection of Inputs, Outputs and a Command that are identified by a _name_. Inputs and Outputs are as above. Note that you can use functions such as `stdout()` or other for the expression.

> If you don't want to play by these rules, don't include any inputs or outputs and just provide your whole string to the initializer for command.

```python
t = wdlgen.Task("task_name")
t.inputs.append(wdlgen.Input(wdlgen.WdlType.parse("String"), "taskGreeting"))
# command in next section
t.outputs.append(wdlgen.Output(wdlgen.WdlType.parse("File"), "standardOut", "stdout()"))
```

#### Command

The command is broken up similar to how CWL breaks its command generation up, by itself it has a _base command_. Each component has a corresponding input (else use the `wdlgen.Task.Command.Argument` class), optionality, position, prefix (and whether the value should be separated from prefix; think `-o {val}` vs `outputDir={val}`) and potentially a default.

Construct a command like the following:
```python
command = wdlgen.Task.Command("echo")
command.inputs.append(wdlgen.Task.Command.CommandInput("taskGreeting", optional=False, position=None, prefix="-a", separate_value_from_prefix=True, default=None))
command.inputs.append(wdlgen.Task.Command.CommandInput("otherInput", optional=True, position=2, prefix="optional-param=", separate_value_from_prefix=False, default=None))

# t is the task
t.command = command
print(command.get_string())
```

This will result in the following WDL command:
```bash
echo \
  -a ${taskGreeting} \
  ${"optional-param=" + otherInput}
```

#### Task output:

The combination of the task and command outputs:
```wdl
task task_name {
  String taskGreeting
  command {
    echo \
      -a ${taskGreeting} \
      ${"optional-param=" + otherInput}
  }

  output {
    File standardOut = stdout()
  }
}
```

### Workflow

You should have moderate idea of the structure of WDL as there's no cleverness or abstraction done anywhere. Beware: there's also no checking attributes (to see if your `inputMap` actually corresponds to inputs).

The structure of a workflow is m

```python
w = wdlgen.Workflow("workflow_name")

w.imports.append(wdlgen.Workflow.WorkflowImport("tool_file", ""))
w.inputs.append(wdlgen.WdlType.parse("String"), wdlgen.Input("inputGreeting"))


inputs_map = {"taskGreeting": "inputGreeting"}
w.calls.append(wdlgen.WorkflowCall("Q.namspaced_task_identifier", "task_alias", inputs_map))
w.outputs.append(wdlgen.Output(wdlgen.WdlType.parse("File"), "standardOut", "task_alias.standardOut")
```

Which outputs:
```wdl
import "tools/tool_file.wdl"

workflow workflow_name {
  String inputGreeting
  call Q.namspaced_task_identifier as task_alias {
    input:
      taskGreeting=inputGreeting
  }
  output {
    File standardOut = task_alias.standardOut
  }
}
```


## Known limitations

I'm not a fan of the string interpolation generation of WDL that this module does. I think trying to build an [Abstract syntax tree](https://en.wikipedia.org/wiki/Abstract_syntax_tree) and then there should be something that convert that into the DSL that WDL uses.

You could also cause syntax errors in generated WDL by providing illegal characters. 

## Goals

- Improve code-level documentation.
- Increase the testing coverage + quality of unit tests.
- Better represent the WDL spec.
- Find an easier distribution / release method - such as PIP.
- Automate testing and delivery through TravisCI / CircleCI or similar.
- Validate each value by [WDL's language specifications](https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md#language-specification).

### Long goals
- Write a documentation site.
- Make classes convert into AST and then into DSL.

## Issues and Pull Requests
Feel free to log issues and make pull requests. I make no guarantee to the existence or timeliness of replies.


Links:

- WDL description: https://github.com/openwdl/wdl/blob/master/versions/1.0/SPEC.md