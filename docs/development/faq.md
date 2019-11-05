# Frequently asked questions

This document contains some common FAQ that may not
have a place outside this documentation.


## Command tools and the command line

- **Why is my input not being bound onto the command line?**

    You need to provide a ``position`` or ``prefix`` to a :class:`janis.ToolInput` to be bound on the command line.

- **How do I prefix each individual element in an array for my command line?**

    Set ``prefix_applies_to_all_elements=True`` on the :class:`janis.ToolInput`.

- **How do I make sure my file is in the execution directory? / How do I localise my file?**

    Set ``prefix_applies_to_all_elements=True`` on the :class:`janis.ToolInput`.
    
- **How do I include environment variables within my execution environment?**

    These are only available from within a CommandTool, and available by overriding the ``env_vars(self)`` method to return a dictionary of string to ``Union[str, Selector]`` key-value pairs.

- **Can a `janis.ToolInput` be marked as streamable?**

    Currently there's no support to mark ToolInput's as `streamable`. Although the
    [`CWL specification`](https://www.commonwl.org/v1.1/CommandLineTool.html#CommandInputParameter).
    does support marking inputs as streamable, WDL does not and, there is 
    [no engine support](https://github.com/broadinstitute/cromwell/issues/3454#issuecomment-455367417). 
      
    