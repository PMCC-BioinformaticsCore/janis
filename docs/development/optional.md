# Optional types

Just a couple of notes about optional types and applied defaults, because CWL and WDl 
handle these differently internally. From Janis they look like they are handled the same.

In Janis, if a type has a default, it is considered _Optional_. If the value of the input 
doesn't resolve to a value, the default is applied. 

For the purposes of this discussion, we're going to borrow two concepts from JavaScript:

- `undefined` - This is where there was no value applied. It means that no value
was ever set and is left in an uninitialised state.

- `null` - This is an assignment value, which can have meaning assigned to it. It often
has the meaning that it's unknown.


## CWL

In CWL, a default value is applied if the value [resolves to null](https://www.commonwl.org/v1.1/Workflow.html#WorkflowStepInput).

> The default value for this parameter to use if either there is no source field, 
or the value produced by the source is null. The default must be applied prior 
to scattering or evaluating valueFrom.

This means that all the way up to the almost end can be a _nullable_ value.

CWL hence makes no distinction between a `null` or `undefined` state. An input with a default
may be _assigned_ an `undefined` value (not mapped) or `null` and the default will be applied. 

## WDL 

> In our discussion, we'll use `null` to reference the [_optional literal_](https://github.com/openwdl/wdl/blob/master/versions/development/SPEC.md#optional-literals)
which in WDL's case is `None`.

In WDL, a default is only applied when the value is not passed, that is the variable is in
an `undefined` state. This means that a variable with a default is not necessarily optional
and may not have a `null` values assigned to it.


- `String q1` - The variable `q1` must be mapped AND have a value assigned to it at runtime.

- `String? q2` - The variable `q2` may have a value or `None` assigned, or simply not be mapped. 

- `String q3 = "test1"` - A variable `q3` may remain unmapped or must be assigned a string. If
there is no mapping, a value of `"test"` will be applied.

- `String? q4 = "test2"` - The variable `q4` may remain unmapped in which case the value `"test2"`
is assigned to it. It may also be passed a value or `None` which will override the default.


Janis uses the [`select_first`](https://github.com/openwdl/wdl/blob/master/versions/development/SPEC.md#x-select_firstarrayx)
or [`defined`](https://github.com/openwdl/wdl/blob/master/versions/development/SPEC.md#boolean-definedx) & 
[`If then else`](https://github.com/openwdl/wdl/blob/master/versions/development/SPEC.md#if-then-else)
operators to provide default values to a variable or argument.

Sometimes, a [`default` expression placeholder option](https://github.com/openwdl/wdl/blob/master/versions/development/SPEC.md#default)
may be used to provide a default command line binding.
 