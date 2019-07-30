# Examples

In this section, we give a few basic examples of how to build
simple workflows with Janis.

> See the [Janis example pipelines repository](https://github.com/PMCC-BioinformaticsCore/janis-examplepipelines)
for whole genome analysis bioinformatics workflows.

## Workflows

### Echo

> Echo a string to the console.

Link: https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/echo.py

This is the same workflow as on the README.

### Simple Workflow

> Untar an archive, compile the contents and re-tar the result.

- [Tutorial guide](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/simple.md):
 
    - [Progressive workflow](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/simple.py):
    The simple workflow progressively constructed in a standard Python file.

    - [Subclass workflow](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/simplewrapped.py): 
    The same workflow, but subclassed to demonstrate how this might be achieved to be easily reused.  

## Functionality Demonstrations

### Edges

> A simple _non-functional_ workflow that demonstrates different ways to connect edges.

https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/edges.py

### Secondary Files 

> A simple _non-functional_ workflow that demonstrates how to define a data type 
with secondary files and how to use it.

https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/secondaryfiles.py

It's easy to see how the automatic generation of secondary inputs happens in the `WDL`
output of this workflow.

