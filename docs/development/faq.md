# Frequently asked questions

This document contains some common FAQ that may not
have a place outside this documentation.

- **Can a `janis.ToolInput` be marked as streamable?**

    Currently there's no support to mark ToolInput's as `streamable`. Although the
    [`CWL specification`](https://www.commonwl.org/v1.1/CommandLineTool.html#CommandInputParameter).
    does support marking inputs as streamable, WDL does not and, there is 
    [no engine support](https://github.com/broadinstitute/cromwell/issues/3454#issuecomment-455367417). 
      
    