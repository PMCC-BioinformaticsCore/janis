# Contributing

We'd love to accept community contributions, especially to add new tools to:

- Janis-bioinformatics
- Janis-unix
- Or other tool suites!

There are a couple of best standard we'd love to follow when contributing to Janis.

## Documentation

There are a couple of ways that can you provide documentation to your tools:

### Tool and workflow inputs

There's a `doc` parameter that can be filled when creating a `ToolInput` or `Workflow.input` that is used to drive the Janis documentation and table of inputs for a tool.

```python
ToolInput(
	"inputId", 
	DataType, 
	doc="<Explain what this input does>", 
	**kwargs
)
```

```python
w = WorkflowBuilder()
w.input(
	"inputId", 
	DataType, 
	doc="<Explain what this input does>", 
	**kwargs
)
```

### Tool documentation

This differs if you're using a class based tool, or using of the builders.

#### Class based

A `CommandTool` and `Workflow` both have a `bind_metadata` method that you can return a `ToolMetadata` or `WorkflowMetadata` from respectively.

**CommandTool**:
```python
def bind_metadata(self):  
	from datetime import date  

	return ToolMetadata(
		contributors=["<Add your name here>"],
		dateCreated=date(2018, 12, 24),
		dateUpdated=date(2020),
		institution="<who produces the tool>",
		keywords=["list", "of", "keywords"],
		documentationUrl="<url to original documentation>",
		short_documentation="Short subtitle of tool",
		documentation="""Extensive tool documentation here """.strip(),
		doi=None,
		citation=None,
  )
```

**Workflow**

```python
def bind_metadata(self):  
	from datetime import date  

	return ToolMetadata(
		contributors=["<Add your name here>"],
		dateCreated=date(2018, 12, 24),
		dateUpdated=date(2020),
		institution="<who produces the tool>",
		keywords=["list", "of", "keywords"],
		documentationUrl="<url to original documentation>",
		short_documentation="Short subtitle of tool",
		documentation="""Extensive tool documentation here """.strip(),
		doi=None,
		citation=None,
  ).   
```

#### Builders

The builders have a parameter called `metadata` which takes a `ToolMetadata` or `WorkflowMetadata` for a `CommandToolBuilder` or `WorkflowBuilder` respectively:

**Example**
```python
wf = WorkflowBuilder(
    "workflowId",
    **kwargs,
    metadata=WorkflowMetadata(
        **meta_kwargs
    ),
)

# OR

clt = CommandTooBuilder(
    "commandtoolId",
    **kwargs,
    metadata=ToolMetadata(
        **meta_kwargs
    ),
)
``` 

## Code formatting

Janis uses the opinionated [Black](https://github.com/psf/black) Python code formatter. It ensures that there are no personal opinions about how to format code, and anecdotally saves a lot of mental power from not thinking about this.

The projects include a [`pyproject.toml` ](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics/blob/master/pyproject.toml), and it's a git pre-commit dependency.

### Installation

You can install Black (+ pre-commit) inside a `janis-*` directory with:

```bash
pip3 install -r requirements.txt
```

You can manually run black with:

```bash
black {source_file_or_directory}
```

### Editor integration

> See more: [https://github.com/psf/black](https://github.com/psf/black)

We'd recommend configuring Black to run on save within your editor with the following guides:

- [PyCharm / IntelliJ IDEA](https://github.com/psf/black#pycharmintellij-idea)
- [Vim](https://github.com/psf/black#vim)
- [VSCode](https://code.visualstudio.com/docs/python/editing#_formatting) + [Python extension](https://marketplace.visualstudio.com/items?itemName=ms-python.python), simplified instructions:
    > - Install the [extension]
    > - Within VSCode, go to your settings.json
    >   - Command + Shift + P
    >   - (type) Preferences: Open Settings (JSON)
    > - Add the following lines:
    >   ```
    >       "python.formatting.provider": "black",
    >       "editor.formatOnSave": true
    >   ```



