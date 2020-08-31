CONTRIBUTING

Thanks for considering contributing to Janis! Whether by providing feedback, raising issues or contributing code, we thank you for engaging and making Janis a better pipeline framework!! 🙏

This document are a set of guidelines when contributing to Janis and its packages. These are mostly guidelines, not rules. Please use your best judgment and consider others when engaging with Janis, and feel free to propose changes to this document in a pull request.  


## I have a question?

Before raising an issue to ask a question, we recommend you check the following resources (as they'll probably be a bit quicker):

- [Frequently asked questions](https://janis.readthedocs.io/en/latest/references/faq.html)
- [Common Errors](https://janis.readthedocs.io/en/latest/references/errors.html)
- [Janis patterns](https://github.com/PMCC-BioinformaticsCore/janis-patterns)
- [Gitter#janis-pipelines](https://gitter.im/janis-pipelines/community)

## Structure

Janis is split into lots of smaller components, this allows us to have smaller incremental releases to iterate and fix bugs quicker - at the expense of greater complexity.

- [Janis (main)](https://github.com/PMCC-BioinformaticsCore/janis) - Documentation, placeholder and other various components.
- [janis-core](https://github.com/PMCC-BioinformaticsCore/janis-core) - Python library for building and translating workflows
- [janis-assistant](https://github.com/PMCC-BioinformaticsCore/janis-assistant) - a layer for automatically configuring and running Janis workflows on CWLTool / Cromwell.
	- [janis-templates](https://github.com/PMCC-BioinformaticsCore/janis-templates) - contains templates for configuring janis-assistant for a number of environments (eg: peter mac, spartan, wehi, pawsey)
- Libraries (prebuilt tools and workflows)
	- [janis-bioinformatics](https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics) ([docs](https://janis.readthedocs.io/en/latest/tools/bioinformatics/index.html)) - bioinformatics tools
	- [janis-unix](https://github.com/PMCC-BioinformaticsCore/janis-unix) ([docs](https://janis.readthedocs.io/en/latest/tools/unix/index.html)) - unix tools
	- [janis-pipelines](https://github.com/PMCC-BioinformaticsCore/janis-pipelines) ([docs](https://janis.readthedocs.io/en/latest/pipelines/index.html)) - prebuilt bioinformatics workflows (primarily WGS pipelines currently)


## Code quality

Janis uses combinations of code-formatting, unit tests and additional scripts to evalute your contribution. These are often run by our CI (Travis), but we recommend running these checks before you open your PR.

### Formatting

All the code in Janis should be formatted using [Black](https://github.com/psf/black), we recommend configuring Black to _on-save_ within your editor. 

- PyCharm / IntelliJ IDEA
- Vim
- VSCode + Python extension, simplified instructions:
	- Install the [extension]
	- Within VSCode, go to your settings.json
		- Command + Shift + P
		- (type) Preferences: Open Settings (JSON)
		- Add the following lines:
		```
	    "python.formatting.provider": "black",
	    "editor.formatOnSave": true
	    ```

### Unit tests

> _Further information: [Testing](https://janis.readthedocs.io/en/latest/development/testing.html?highlight=nosetests)_

Janis-core has decent unittest coverage, and we endeavour to increase testing across the various modules of Janis. 

If you find a bug, please add a unit test to cover your change to decrease likelihood of a regression. 


### Tool definitions

Tool definitions (wrappers, interfaces) are fairly difficult to comprehensively test, and there's not usually a "correct" way to do it. We trust the tool definitions

We're working on ways to automate evaluting some aspects of a tool definition, but please make sure:

- It's correct,
- It translates to BOTH cwl and wdl,
- It's comprehensive,
- It contains: friendly_name, version, metadata (this is important for the documentation generator).


## Expectations

From us, you can expect feedback on PR and issues in a timely manner. There is a chance we might have missed your issue (sometimes GitHub notifications aren't always set up correctly) - please feel free to reply to an issue (to bump it), or tag (@illusional) for any concerns.

From you, we expect you to follow this contribution guide, and act in good faith when engaging with Janis and other contributors.


## Closing notes

We may update this contributing guide from time-to-time with new resources or standards we expect from contributions. Please make sure you check this document from time-to-time.


### Attribution

This contributing guide was explicitly based on the [Atom](https://github.com/atom/atom/blob/master/CONTRIBUTING.md) contribution guide, and implicitly based on other contributing guides on GitHub. 

Janis as a project is based on a large body of previous work which we can directly reference as modules, or is implicitly used by the ideas that Janis captures. We recognise the intellectual effort that Janis is built derived from.
