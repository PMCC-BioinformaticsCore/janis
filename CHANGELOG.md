# Change-log

## v0.9.0

### Upgrade notes

Please read this (incompatible changes) section carefully:

> v0.9.0 brings the following backwards incompatible changes:
> - Workflows:
>   - `Workflow.output` param `output_tag` renamed to `output_folder`
>   - `Workflow.output` param `output_prefix` renamed to `output_tag`
> - Bioinformatics:
>   - Bam index pattern `^.bai` refactored to `.bai`
>     - ie: `filename.bai` → `filename.bam.bai`
>     - see `secondaries_present_as` to fall back to old method
>     - GATK tools have been fixed to present secondary index as `^.bai`.
>     - Most HTSLib should be fine.
> - Assistant
>   - Require an output folder to be specified
>     - Ability to specify output folder multiple times
>     - Execution directory falls under this 
>   - Respect changes to `output_name` + `output_folder` required minor schema change
>  - Templates
>    - Refactor of all templates to use [camel case](https://en.wikipedia.org/wiki/Camel_case). 

### Janis-core changelog

- New: [**PythonTool**](https://janis.readthedocs.io/en/latest/references/tools/pythontool.html)
	- Allow arbitrary, self-contained Python scripts to be executed inside a container.
	- Includes base for other CodeTools.
- Workflows:
  - `Workflow.output` param `output_tag` renamed to `output_folder`
  - `Workflow.output` param `output_prefix` renamed to `output_tag`
 - New: `presents_as` to allow an input to be localised as a specific filename.
- New `secondaries_present_as`, to allow a secondary / accessory file to be presented as a different pattern
	- Requires [common-workflow-language/cwltool/issues/1232](https://github.com/common-workflow-language/cwltool/issues/1232) to be resolved to work correctly.
- Fix `mv` statements in WDL to use `-n` instead of `-t`
- Add new LogLevel verbose, and change default to DEBUG
- Removes dynamic filename generation preparing for better call-caching
- WDL to use `~{}` instead of `~{}` with `command <<< ~{} >>>` syntax. WDLGen: [`92470d7`](https://github.com/PMCC-BioinformaticsCore/python-wdlgen/commit/92470d7acd0ef870210d4891bfc14c9d570e6f07)
- Add better ability to trace JanisShed / registry

### Janis-assistant changelog

#### Big changes

- Engines by default run in the foreground
- Allow ability to specify output directory multiple times
	- Existing functionality where new output directories are created for each run can be restored by adding `output_dir` into your janis config (default: `~/.janis/janis.conf`)
- `janis watch $wid`:
	- Doesn't clear your current history (opens curses display with `blessed`)
	- Ability to specify the output directory
- Change default execution directory to be inside the output directory
- Rename template names to respect camel case
- CLI changes:
	- `--should-disconnect` renamed to `--background / -B`
	- `--use-mysql` renamed to `--mysql`
	- `--no-watch` is default, `--progress` to show progress by default.
	- Added `--development` flag:
		-  sets `--mysql` and `--keep-intermediate-files` to true
- Allow saving of literal output values on a local filesystem

#### Other bug fixes
- Fully qualify output directory when specified on the CLI.
- Check if container is available before run is issued
- Allow CodeTool to be passed in via CLI
- Ensure `janis init --stdout` actually logs to stdout
- More camel-case refactor
- 


#### Engines

- Engines callback the progress by default, rather than Janis polling at some arbitrary interval,
- Cromwell polling logic exists inside the engine now
- This is to remove workflow job state from engine, but also prepare for other engines that have the better ability to callback.
- Added ability to configure CWLTool
	- Add singularity support for CWLTool for all Singularity based templates
- Add better support for catching Cromwell + CWLTool exit errors
- Cromwell uses the localization strategy `[hard-link, cached-copy]` ([source](https://github.com/broadinstitute/cromwell/pull/4900#issuecomment-571456289)) 

Notes about Cromwell call-caching:
- More effort has been added to support Cromwell's call-caching, however the following issues require resolving:
  - [broadinstitute/cromwell/#5348](https://github.com/broadinstitute/cromwell/issues/5348) [[BA-6172](https://broadworkbench.atlassian.net/browse/BA-6172)]
  - [broadinstitute/cromwell/#5346](https://github.com/broadinstitute/cromwell/issues/5346)


## v0.8.x

### Upgrade notes:

> On a CommandTool, you must:
> - Remove @staticmethod and add (self) as a parameter to the following methods on a tool: tool(), container(), version(), tool_module(), tool_provider() [eg to become container(self)].


### Documentation

- Add documentation for new features such as CommandToolBuilder

### Janis-core changelog

This allows the creation of a `CommandToolBuilder` like the following:

- New CommandToolBuilder for simplified building of command tools:
	- This tool syntax should be interchangeable with the existing class inherit syntax. 

```python
ToolName = j.CommandToolBuilder(
    tool: str="toolname",
    base_command=["base", "command"],
    inputs: List[j.ToolInput]=[],
    outputs: List[j.ToolOutput]=[],
    container="container/name:version",
    version="version",
    friendly_name=None,
    arguments=None,
    env_vars=None,
    tool_module=None,
    tool_provider=None,
    metadata: ToolMetadata=j.ToolMetadata(),
    cpu: Union[int, Callable[[Dict[str, Any]], int]]=None,
    memory: Union[int, Callable[[Dict[str, Any]], int]]=None,
)
```

### janis-assistant changelog
Primarily the new methods are:

- `spider` look up documentation about an existing tool
- `docs` opens the docs in your webbrowser, this works for browsers without an associated web browser

Quality changes to user experience when using Janis assistant:

- Newer + extensible templates


## v0.7.x

### Documentation
- New sections (Tutorials, Registry, References, Development)  
- More complete tutorial for setup + wrapping new tool  
- New guide for containerising tools  
- Basic CWL equivalency guide  
- FAQ: [https://janis.readthedocs.io/en/latest/references/faq.html](https://janis.readthedocs.io/en/latest/references/faq.html)  
- Common errors: [https://janis.readthedocs.io/en/latest/references/errors.html](https://janis.readthedocs.io/en/latest/references/errors.html)  
- Updated simple workflow: [https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/simple.md](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/examples/simple.md) (edited)

### Janis-core changelog

- New: Scattering by multiple fields (dot + cross) for CWL and WDL
- New: extensible input validation with errors reported back to user
- New: stderr type
- Fixes: Logger truncating logfile on reload

### Janis-assistant changelog 

> This is a big and backwards **incompatible release**

Migration guide:

```
# Remove config and task DB
rm ~/.janis/janis.conf
rm ~/.janis/janis.db
janis init local
```

In this exciting (and unfortunately backwards incompatible) release of Janis, Janis now takes care of metadata received by the engines. It doesn't require Cromwell to be open to poll the status of the workflow, though Cromwell metadata is still used to update its own store.

Now that execution and management of a workflow can be run in the background, the workflow "brain" (Janis + Cromwell) can be submitted to the cluster to continuously monitor Cromwell (also running on a cluster), and means that notifications can be sent all the time.

The config files and existing templates have been changed, your existing templates might not work correctly. It's advised to remove your existing config `rm ~/.janis/janis.conf` and recreate it with `janis init local`. 

Other features ([Commit diff](https://github.com/PMCC-BioinformaticsCore/janis-runner/compare/v0.6.2...v0.7.1)):

- New metadata DB store
- Output prefix + tag support through Janis annotations
    - Group your outputs into folders or by prefix with these annotations
- Better watching of Cromwell with disconnected `janis monitor`
- Stability with Slurm `afternotok` job dependency
- Email notifications 
- Better logging to disk of cromwell engine and janis during disconnected execution
- Input validation 
- Help text on CLI
- New templates that can do more
- Better time duration formatter (`n seconds` to `days:minutes:hours:seconds`)

#### Additional release notes:

CLI methods:

-  New  `janis docs`  to take you to the docs - anywhere ;)
-  Revised argument groups for  `janis run -h`
- Changed default Cromwell execution directory to be in `$HOME/execution`  
- Changed default output directory to be `$HOME/janis/`  
- Overridable by three new env vars: `JANIS_BASEDIR`  `JANIS_EXCECUTIONDIR`  `JANIS_OUTPUTDIR`  
- `-o` shortcut for output translation  
- Outputs on the progress screen  
- More stability in CWLTool implementation + abort  
- More descriptive error messages - relayed in docs
Templates can now be extended through the  `janis.templates`  extension. See the new repo:  [https://github.com/PMCC-BioinformaticsCore/janis-templates](https://github.com/PMCC-BioinformaticsCore/janis-templates)  for more information on editing existing and creating your own templates.

Stability changes for MariaDB for singularity ([docker-library/mariadb/issues/262](https://github.com/docker-library/mariadb/issues/262))



## v0.6.0

New Janis API!

> **Warning**: Janis v0.6.0 brings backwards incompatible change to workflow building.

In this exciting release, we present a new revised Janis API for workflow construction. Instead of creating nodes and adding them to a workflow, we can use the new `WorkflowBuilder` to simplify the process of connections directly within the step. Here's an example from the readme:

```python
import janis as j
# Nb: in v0.7.0 this becomes `from janix.tools import Echo`
from janis.unix.tools.echo import Echo

w = j.WorkflowBuilder("workflowId")

w.input("inputIdentifier", j.String, default="Hello, World!")
w.step("stepIdentifier", Echo(inp=w.inputIdentifier))
w.output("outputIdentifier", source=w.stepIdentifier.out)
```

There are new tools, updated docs and more examples within janis to keep you moving quicker. We'd like to invite you to join the Janis Gitter: https://gitter.im/janis-pipelines/community where you can ask any questions about Janis.

There are also a number of improvements to running workflows using `janis-runner`, now you can run a workflow as simple as:

```bash
janis run hello --engine [cwltool|cromwell] --inp "Hello, Janis!"
```



## v0.5.x

### v0.5.5


### v0.5.1

Commit: [`0bda65b`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.5.0)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis-core/compare/v0.5.0...v0.5.1)

- Better inheritance checks for File subclasses
- Cleaner`translate` API for `CommandTool`
- Autowrap tool in workflow `CommandTool.wrapped_in_wf(self)`
- Better documentation on a CommandTool
- Fixed input overrides for CWL files and directories

### v0.5.0

Commit: [`0bda65b`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.5.0)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis-core/compare/v0.4.1...v0.5.0)

- CommandTool API change:
	- Rename `docker() -> container()`
- Support reverse add for `InputSelectors`

## v0.4.x

### v0.4.1

Commit: [`afe7aef`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.4.1)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis-core/compare/v0.4.0...v0.4.1)

-   Ensure defaults for inputs are propagated through InputSelectors within a StringFormatter in WDL.
-   Ensure imports in WDL are not appearing multiple times.
-   Don't create an output folder when  `to_disk`  is False (this changes the behaviour that job files are not generated when  `to_disk=False`.
-   More consistent debugging messages when WDL input selectors are incorrect.

### v0.4.0

Commit: [`afe7aef`](https://github.com/PMCC-BioinformaticsCore/janis-core/releases/tag/v0.4.1)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis-core/compare/v0.3.1...v0.4.0)

There's a new modular structure to Janis! The core functionality has been moved to  `janis-core`, and the unix tools also being split out into  `janis-unix`. These are still dot accessible through  `janis.bioinformatics`  or  `janis.unix`  through use of the new  `janis.extension`  [entry point](https://amir.rachum.com/blog/2017/07/28/python-entry-points/). In theory, this allows anyone to add their own set of tools and extensions.

There's no resulting change for the user, all of the tests pass with no changes, so please raise an issue if you're having troubles.

This release sees tools REQUIRING a version tag now, and the docs registering the use of different versions of tools and documentation (there's also an updated style guide).

> All tools are synchronised to v0.4.0 or greater.

## v0.3.x

### v0.3.3

Commit: [`805e6a3`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.3.3)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/805e6a3...c1a0943)

-   New:  `extension`  property on the  `File`  data type for the  _expected_  extension as a hint for the secondary file collection.
-   Quality: Ability to reverse add StringFormatters to a string, ie:  `"test_" + StringFormatter(*kwargs)`(coming soon, reverse add for other selectors.
-   Cascading defaults in WDL (to reflect similar behaviour in CWL) by internally utilising the WDL  `select_first`  and  `default`  functions. This requires an update to  `illusional.wdlgen==0.2.3`.
-   CPU now has a non-optional default of 1 in WDL.

### v0.3.2

Commit: [`805e6a3`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.3.2)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/5aaeb2d...805e6a3)

- New: `extension` property on the `File` data type for the _expected_ extension as a hint for the secondary file collection.
- Quality: Ability to reverse add StringFormatters to a string, ie: `"test_" + StringFormatter(*kwargs)` (coming soon, reverse add for other selectors.
- Cascading defaults in WDL (to reflect similar behaviour in CWL) by internally utilising the WDL `select_first` and `default` functions. This requires an update to `illusional.wdlgen==0.2.3`.
- CPU now has a non-optional default of 1 in WDL.


### v0.3.1

Commit: [`5aaeb2d`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.3.1)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/952e85c...5aaeb2d)

- Changes export location to current working directory
	- Adds more stability for exports
- Converts some connection warnings into "info" statements
- Fully qualifies all connections in the examples

### v0.3.0

Commit: [`952e85c`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.3.0)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/dfc13bc...952e85c)

This commit brings some backwards incompatible changes. 

New:
- String Formatters
  - Basic syntax: `StringFormatter("format including {arg}", arg=InputSelector("valueToSelect"))`
  - String formatters, InputSelectors and Strings can be concatenated to create more complex string formats. (Nb: If you want a prefix, you'll need to wrap a string in an empty string formatter: ie:
    - Ok:
      - `StringFormatter("value") + "suffix"`
      - `StringFormatter("prefix") + StringFormatter("value")`
    - Bad:
      - ~~`"prefix" + StringFormatter("value")`~~
  - For this reason, the prefix and suffix components have been removed from InputSelector.
- Generate a table of resources a workflow


Improved:

- Upgraded CWL support (working with CWLTool)
  - Automatic ShellCommandRequirement
  - Support for prefix on every element of optional array.
- Better errors 
- More tests across translations and functionality

## v0.2.x

### v0.2.17
Commit: [`dfc13bc`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.2.17`)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/8050495...master)

- Add [Command Line Tool generator](https://github.com/PMCC-BioinformaticsCore/janis/blob/master/janis/utils/parse_clt.py) to assist in the tool-wrapper generation process.
- Improve tests and increase test coverage

	- This avoids the accidental overriding when generating inputs (to avoid having to remove default values).
	- This attribute can only apply to optional or defaulted inputs.
- Adds JSON and CSV File types
- Fix WDL workflow validation, looks for `$womtool` in validation.
- Ensures stepIds are case sensitive
- Begin the streamlining process for [Shepherd](https://github.com/PMCC-BioinformaticsCore/shepherd)
	- Adds `include_in_inputs_file=True` param to `Input`, to avoid including an input in a file if the provided value is `None`. 	
	- Value checking on translation (`allow_null_if_not_optional`)
- More enums instead of _arbitrary_ string values


### v0.2.16
Commit: [`8050495`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.2.16`)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/4057635...8050495)

- Ensure `runtime_disk` is not generated at a tool level when resource overrides are not generated.
	- `runtime_memory` and `runtime_cpu` are always generated as some tools may rely on these values through the CPU and Memory selectors.

### v0.2.15
Commit: [`4057635`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.2.15`)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/bc7f8ad...4057635)

- Increases testing coverage
- Better example codes within docs and example folder.
- Memory and CPU selector defaults
- Removes adding defaults to edges:
	- The preferred method is adding a default to an `Input`, this ensures cascading defaults
	- Note, passing a nulled value to an input will override, even if there's a default.
	- This highlights the the runtime_disk attribute which is not nullable, even though the value is optional.

### v0.2.14
Commit: [`bc7f8ad`](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.2.14`)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/732d5dd...bc7f8ad)

- Removes default from DataType block in favour of on `Input` and `ToolInput`
- Improve code quality
- Add CWL Resource overrides

### v0.2.13

Commit: [732d5dd](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.2.13)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/10b69ad...732d5dd)

- Reopens library with new GPL3 License. 
-   More helpful error messages
-   More and better support for Selectors,
    -   Add more mechanisms for referencing resources within your tools: CpuSelector and MemorySelector (only WDL at the moment)
-   Stronger foundation for WDL resource overrides
-   Better localisation of files across CWL and WDL

### v0.2.12

Commit: [10b69ad](https://github.com/PMCC-BioinformaticsCore/janis/releases/tag/v0.2.12)
Code level: [Comparison](https://github.com/PMCC-BioinformaticsCore/janis/compare/4d77eee...10b69ad)

- More functional CWL and WDL
- Better secondary file support with scatters in WDL through the use of `transform` .
	- Ensure that you're only specifying the extension on generated filenames, and using the `suffix` field.

## Other releases

Releases before here aren't tracked as well, but they are still present on the Github [Janis/releases](https://github.com/PMCC-BioinformaticsCore/janis/releases) page.