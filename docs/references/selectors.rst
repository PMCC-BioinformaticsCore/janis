Selectors and StringFormatting
==============================

*Methods to reference or transform inputs*

Overview
*********

A :class:`janis.Selector` is an abstract class that allows a workflow author to
reference or transform a value. There are 3 main forms of selectors:

- :class:`janis.InputSelector` - Used for selecting inputs
- :class:`janis.WildcardSelector` - Used for collecting outputs
- :class:`janis.StringFormatter` - Constructing or transforming strings


InputSelector
..............

Declaration
-----------

.. autoclass:: janis.InputSelector

Overview
--------

An ``InputSelector`` is used to get the value of an input at runtime.
This has a number of use cases:

- Collecting outputs through the ``glob=`` field.
- Constructing a :class:`janis.StringFormatter` from existing inputs

These inputs are checked for validity at runtime. In CWL, defaults are propagated by
definition in the language spec, in WDL this propagation is handled by janis.

Example
-------

The following code will use the value of the ``outputFilename`` to find the corresponding
output to ``bamOutput``. Additionally as the data_type ``BamBai`` has a secondary file,
in WDL janis will automatically construct a second output corresponding to the ``.bai``
secondary (called ``bamOutput_bai`` and modify the contents of the ``InputSelector``
to correctly glob this index file.

.. code-block:: python

   def outputs(self):
       return [
           ToolOutput("bamOutput", BamBai, glob=InputSelector("outputFilename"))
       ]


WildcardSelector
.................

Declaration
-----------

.. autoclass:: janis.WildcardSelector

Overview
--------

A wildcard selector is used to glob a collection of files for an output where
an ``InputSelector`` is not appropriate. Note that if an ``*`` (asterisk | star)
bind is used to collect files, the output type should be an array.

Example
-------

The following example will demonstrate how to get all the
files ending on ``.csv`` within the execution directory.

.. code-block:: python

   def outputs(self):
       return [
           ToolOutput("metrics", Array(Csv), glob=WildcardSelector("*.csv"))
       ]


StringFormatting
.................

Declaration
-----------

.. autoclass:: janis.StringFormatter

Overview
--------

A ``StringFormatter`` is used to allow inputs or other values to be inserted at runtime
into a string template. A ``StringFormatter`` can be concatenated with Python strings,
another ``StringFormatter`` or an ``InputSelector``.

The string ``"{placeholdername}"`` can be used within a string format, where ``placeholdername``
is a ``kwarg`` passed to the StringFormatter with the intended selector or value.

The placeholder names must be valid Python variable names (as they're passed as kwargs).
See the `String formatter tests <https://github.com/PMCC-BioinformaticsCore/janis-core/blob/master/janis_core/tests/test_stringbuilder.py>`_
for more examples.

Example
--------

- Initialisation

    ``StringFormatter("Hello, {name}", name=InputSelector("username"))``

- Concatenation

    ``"Hello, " + InputSelector("username")``
    ``InputSelector("greeting") + StringFormatter(", {name}", name=InputSelector("username"))``
