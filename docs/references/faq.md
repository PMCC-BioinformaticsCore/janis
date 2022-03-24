# Frequently asked questions

This document contains some common FAQ that may not
have a place outside this documentation.


## Janis Registry

- **Exception: Couldn't find tool: 'MyTool'**

To ensure your tool is correctly found by the registry (`JanisShed`), you must ensure that:

- The tool implements all abstract methods (can be initialised)
- Available within two levels of the janis-pipelines.tool extension root. This means that it needs to be within one additional import from the base `__init__`.

For example:

```
.
|-- __init__.py     # imports my_directory
|-- my_directory/
|   |-- __init__.py # imports tool from mytool.py
|   |-- mytool.py
```

If you're still having trouble, use `janis spider --trace mytool` to give you an indication of why your tool might be missing.


## Command tools and the command line

- **Why is my input not being bound onto the command line?**

    You need to provide a ``position`` or ``prefix`` to a :class:`janis.ToolInput` to be bound on the command line.

- **How do I prefix each individual element in an array for my command line?**

    Set ``prefix_applies_to_all_elements=True`` on the :class:`janis.ToolInput`.

- **How do I make sure my file is in the execution directory? / How do I localise my file?**

    Set ``localise_file=True`` on the :class:`janis.ToolInput`. Although the ``localise_file`` attribute is allowed for Array data types, the WDL translation will become disabled as this behaviour is not well defined.
    
- **How do I include environment variables within my execution environment?**

    These are only available from within a CommandTool, and available by overriding the ``env_vars(self)`` method to return a dictionary of string to ``Union[str, Selector]`` key-value pairs.

- **Can a `janis.ToolInput` be marked as streamable?**

    Currently there's no support to mark ToolInput's as `streamable`. Although the
    [`CWL specification`](https://www.commonwl.org/v1.2/CommandLineTool.html#CommandInputParameter).
    does support marking inputs as streamable, WDL does not and, there is 
    [no engine support](https://github.com/broadinstitute/cromwell/issues/3454#issuecomment-455367417). 
      
- **How can I manipulate an input value?**

    Janis only provides a limited way of manipulating input values, internally using a [`StringFormatter`](https://janis.readthedocs.io/en/latest/references/selectors.html#stringformatting).
    
    Eg:
    
    - Initialisation
        - ``StringFormatter("Hello, {name}", name=InputSelector("username"))``

    - Concatenation

        - ``"Hello, " + InputSelector("username")``
        - ``InputSelector("greeting") + StringFormatter(", {name}", name=InputSelector("username"))``      
      
      
- **How do I specify an input multiple times on the command line?**

    - Create a `ToolInput` for the input you want to create.
        - Omit any binding options, ie NO `position` and NO `prefix`. 
    - Create a `ToolArgument` with the value `InputSelector("nameofyourinput") and the binding options:
    
    For example:
    
    ```python
    class MyTool(CommandTool):
        def inputs(self):
            return [ToolInput("myInput", String)]
      
        def arguments(selfs):
            return [
                ToolArgument(InputSelector("myInput"), position=1),
                ToolArgument("transformed-" + InputSelector("myInput"), position=2, prefix="--name")
            ]
    ```

      
- **How do I ensure environment variables are set within my execution environment?**

    You can include the following block within your CommandTool:

    ```python
    # Within CommandTool
  
    def env_vars(self):
       return {
           "myvar1": InputSelector("myInput"),
           "myvar2": "constantvalue"
       }
    ```
   
- **How do I make my generated filenames unique for a scatter or from different runs betweens tools?**

    Long story short, you can't. 
    
    But there are a fun few reasons why that's currently the case:
    
    - The filenames are generated at transpile time for a tool wrapper. This means that tasks that use this tool will get the same filename, including if your workflow scatters over this task.
    - WDL doesn't really have a mechanism for achieving dynamic or generated components like this, and CWL was flaky at best.
    - Call caching in Cromwell relies on the command line being constructed, so the generated filenames currently break this, and a purely randomly generated filename would breat this further.
    - We have cascaded filename components on the roadmap, so your filename can be built from a collection of inputs (sort of possible anyway).


## Running workflows

- **How do I detect a workflow's status from the command line?**

    You can request metadata for a workflow You can use the following bash command to detect the status. 

    ```bash
    janis metadata <dir / wid> |grep status|tr -s ' '|cut -d ' ' -f 2
    ```

    Note that if Janis is killed suddenly, it might not have enough time to mark the workflow as failed or terminated. It's worth checking the last_updated_time, if a workflow hasn't been updated in a couple of hours, it's likely Janis has been killed.

    The metadata will have the following keys:

    - https://github.com/PMCC-BioinformaticsCore/janis-assistant/blob/master/janis_assistant/data/enums/workflowmetadatakeys.py

    And the workflow status can have the following states:

    - https://github.com/PMCC-BioinformaticsCore/janis-assistant/blob/master/janis_assistant/data/enums/taskstatus.py

- **How can I override a specific tools' resources, like memory / ram, CPUs, or runtime (coming soon)?**

    Long story short, lookup the override key with:
    
    ```bash
    janis inputs --resources <workflow>  
    ```
  
    It looks like `yourworkflow_embeddedtool_runtime_memory: 8`. More information:  [Configuring resources (CPU / Memory)](https://janis.readthedocs.io/en/latest/references/resources.html).
    
## Containers

- **How do I override the container being used?**

    You can override the container being used by specifying a `container-override`, there are two ways to do this:
    
    - CLI with the syntax: `janis [translate|run] --container-override 'MyTool=myorganisation/mytool:0.10.0' mytool`
    - API: include container override dictionary: `mytool.translate("wdl", container_override={"mytool": "myorganisation/mytool:0.10.0"})`, the dictionary has structure: `{toolId: container}`. The toolId can be an asterisk `*` to override all tools' containers.
