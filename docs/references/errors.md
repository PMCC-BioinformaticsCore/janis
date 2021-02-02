# Common errors

From time to time you'll come across error messages in Janis, and hopefully they give you a good indication of what's happened, but here we'll address some.

If you've got a strange error message and don't know what's happening, post a message on [Gitter]() or raise an issue on GitHub!


## Command line positions

Command line positioning can be a little confusing. When running or translating you should follow this format:

```bash
janis run [run-arguments] yourworkflow.py [workflow-inputs]
```

Where 
- `run-arguments` are options like `--inputs`, `--recipe`, `--hint`, `--engine`, `--stay-connected`, etc.

- `workflow-inputs` are additional inputs you 

Errors:

- **There were unrecognised inputs provided to the tool "yourtool", keys: engine, recipe, hint, otherkey**

    1. In the command line, you might have accidentally placed a run argument in the workflow inputs section. Eg: :
       - Wrong: `janis run yourtool --engine cwltool`
       - Correct: `janis run --engine cwltool yourtool`

    2. You might have also included workflow inputs that yourtool did not recognise. Please ensure that yourtool accepts all the unrecognised keys.

## Finding files

Janis looks in a number of places outside just your current directory to find certain types of files. Janis will look in the following places for your file:

- The full path of the file (if relevant)
- Current directory
- `$JANIS_SEARCHPATH` (if defined in the environment).

If your path starts with `https://` or `http//`, the file will be downloaded to a cache in the configuration directory (default: `~/.janis/cache`) and will not go through the search path procedure.

You can run your command with `janis -d yourcommand` and look for the following lines which tell you how Janis is searching for files:

```
[DEBUG]: Searching for a file called 'hellop'
[DEBUG]: Searching for file 'hellop' in the cwd, '/path/to/cwd/'
[DEBUG]: Attempting to get search path $JANIS_SEARCHPATH from environment variables
[DEBUG]: Got value for env JANIS_SEARCHPATH '/path/to/janis/', searching for file 'hellop' here.
[DEBUG]: Couldn't find a file with filename 'hellop' in any of the following: full path, current working directory (/path/to/cwd/) or the search path.
[DEBUG]: Setting CONSOLE_LEVEL to None while traversing modules
[DEBUG]: Restoring CONSOLE_LEVEL to DEBUG now that Janis shed has been hydrated
# It's now searching the registry and will fail after it can't find it.
```

Errors:

- **Exception: Couldn't find workflow with name: myworkflow**

	This might occur during a `run` or `translate`. Janis couldn't find `myworkflow` in any of the usual spots, nor in the tool registry.

    - If you're referencing a Python file (eg: `myworkflow.py`), try giving the full path to the file or ensure that it's correctly in your search paths.
    - If you're referencing a tool that should be in the store, check the spelling of the tool. The name you reference needs to be the ToolId (under `def tool(): return "toolid"`).
    - If you're referencing a tool that you're putting into the store, ensure there are no syntax / runtime errors in the Python file. You can double check this a few ways:
        - Running `janis translate /path/to/your/file.py wdl` and ensure that runs correctly.
        - Adding an `if __name__ == "__main__": YourWorkflowClass().translate("wdl")` to the bottom of your python file, then running `python /path/to/your/file.py` and checking that it translates correctly.


- **FileNotFoundError: Couldn't find inputs file: myinput.yml**
    
    This is likely because your inputs file `myinput.yml` couldn't be found in the search path.  Try giving the full path to your file.

- **Unrecognised python file when getting workflow / command tool: hello.py :: $PYTHONERROR**

    This error results in a couple of ways:
    
    - There was a problem when parsing your file (`hello.py`). The error (given by `$PYTHONERROR`) was returned by Python when trying to interpret your file and get the workflow out. This could be syntactic or runtime error in your file.
    
    - You were looking a tool up in the registry, but you had a file / folder in your search path that caused a clash. You can include the `--registry` parameter after the `run` to only look up the tool in the registry.

- **There was more than one workflow (3) detected in 'hel.py' (please specify the workflow to use via the `--name` parameter, this name must be the name of the variable or the class name and not the workflowId). Detected tokens: 'HelloWorkflow' (HelloWorkflow), 'HelloWorkflow2' (HelloWorkflow2), 'w' (HelloWorkflow)**

    When looking for workflows and command tools within the file `myworkflow.py`, Janis found more than one workflow / command tool that could be run, in this case there were 3:
    
    - *HelloWorkflow* (`HelloWorkflow`)
    - *HelloWorkflow2* (`HelloWorkflow2`)
    - *w* (`HelloWorkflow`)

    You need to specify the `--name` parameter with one of these ids (`HelloWorkflow`, `HelloWorkflow2` or `w`) to run / translate / etc.
    
- **ValueError: There were errors in 2 inputs: {'field1': 'value was null', 'field2': 'other reason'}**

    One or more of your inputs when running were invalid. The dictionary gives the field that was invalid (key) and the reason Janis thinks your value was invalid (value).


## Containers

- **Exception: The tool 'MyTool' did not have a container. Although not recommended, Janis can export empty docker containers with the parameter `allow_empty_container=True` or `--allow-empty-container`**

    One of the tools that you are trying to run or translate did not have a container attached to this. In Janis, we strongly encourage the use of containers however it's not always possible. 
    
    To allow a tool to be translated correctly, there are two mechanisms for doing so:
    
    - Translate without a container
    
        - CLI: Include the `--allow-empty-container` flag, eg: `janis [translate|run] --allow-empty-container mytool`
        - API: Include `allow_empty_container=True`, `mytool.translate("cwl", allow_empty_container=True)`.
        
    - Override the container
    
        - CLI with the syntax: `janis [translate|run] --container-override 'MyTool=myorganisation/mytool:0.10.0' mytool`
        - API: include container override dictionary: `mytool.translate("wdl", container_override={"mytool": "myorganisation/mytool:0.10.0"})`, the dictionary has structure: `{toolId: container}`. The toolId can be an asterisk `*` to override all tools' containers. 