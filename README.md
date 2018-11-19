
# WEHI pipeline definition language
A framework for creating specialised, simple workflow definitions that are then converted to Common Workflow Language definitions.

## Python dependencies:

Python dependencies are listed in requirements.txt

Local environment can be initialized by the following command:

>pip install -r requirements.txt

To update the requirements.txt to reflect the latest dependency requirements,

>pip freeze > requirements.txt


Uses NetworkX library: https://networkx.github.io/

## Requirements

This document should consume the WEHI Pipeline Definition Language and emit a translation into CWL (later: WDL). It should contain:
    - A fully specified workflow.
    - A zip of the references (the tools).
    - An input yaml file.

It should have the option to produce a visual graph (as it build the DAG) and perform basic typechecking between connections.

## WEHI Pipeline Definition Language
Following discussions with Evan, I've put together a little _guide_ on how I've interpreted the pipeline language (_name?_).

Theoretically we should be able to consume any input format that will serialize as a dictionary (YAML, JSON), but we'll stick to YAML for the descriptions.

**NB:** There is no more type inferencing if you're seeing this message.

### Linking inputs
You must specify and exactly link inputs to the pipelines, and will likely reference a step's output in another input, you should know what files your tool exports.


```yaml
inputs:                                 # Dictionary of input types, input_labels must be unique
    $input_label:
        $input_type:
            property1: value
            property2: value

    $input_label2:
        $input_type2:
            property1: value
            property2: value

outputs:
    # TBA

steps:    # dictionary of steps
    $step_label:
        tool: $tool/version             # Should be able to just say 'tool' or 'toolCategory' as well
        inputs:
            $input1: $input_label
            $input2: $input_label2

    $step_label2:
        tool: $tool_category
        $input1: $step_label/output1    # You must refer to a tool's documentation to find out the types it exports
        $input2: $input_label2
```

For a very simplified step, I'd be happy to just specify the tool, and all outputs from the previous step is provided to the current step, ie. a purely linear workflow, eg:
```yaml
steps:
    $step1: $tool1
    $step2: $tool2
    $step3: $tool3
```

Notes:
- You must NOT have a label called "input"

### Questions

- Inheritance with type checking
- How are secondary files dealt with, are they all passed?