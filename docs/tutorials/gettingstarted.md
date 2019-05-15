
# Getting started

This guide will get you started using `janis`. 

### Very quick start:
You've already got Python (>= 3.6) installed, just run:
```bash
pip3 install janis-pipelines
```
Import `janis` in python as:
```python
import janis
```

We'd recommend you use a Python IDE like PyCharm for code completion, syntax help and auto-importing.


## A bit longer

Okay, let's walk through things at a more reasonable pace. Let's talk about the individual components that `janis` requires to run on your computer.

### Unix computer

Workflow tools (especially containers) are often only developed for nix, or that's where they run smoothest. Although a lot of the software likely works on Windows, we'd recommend a unix system (that includes macOS) for the best experience.

### Python

Python is a powerful programming language that's very common within the scientific community. Janis is developed in Python and it's also how you can construct workflows, so you need at least Python 3.6.

> Although Janis requires Python to build workflows, once you have a workflow description language, you can run your workflows anywhere your workflow engine will run. You can find more information about that on the [running workflows](/tutorials/runningworkflows) tutorial.

### janis

You can easily install `janis` through PIP (and PyPi) by running:
```bash
pip3 install janis-pipelines
```

If you want to install the bioinformatics tool suite as well, you should run:
```bash
pip3 install janis-pipelines[bioinformatics]
```

[//]: # (Conda has been proposed, and will be available after the repo has been configured)
[//]: # (```bash)
[//]: # (conda install janis)
[//]: # (```)

## Importing into Python

You can import `janis` into python with the following statement:

```python
import janis
```

If you only require a subset of `janis` subclasses, you can use the following statement, followed by any other classes you want to import:

```python
from janis import Workflow, Input, Step, Output
```

## Additional recommendations

### Python IDE
We'd recommend you use a Python IDE like PyCharm as it gives you some awesome features such as code completion, immediate syntax help and auto-importing (this is super useful for the bioinformatics tedious paths).

### Docker
You can build build Workflows with just the bare essentials, but to run you'll need to use a container engine. For local environments (if you're allowed), we'd recommend you install Docker to run the tools. Using containers when you run your workflows allows you to be confident that your tools will act in predictable and _reproducible_ ways when you run them in other environments.

## Finished!

Congratulations on getting `janis` installed, and imported into Python! You can now move on to writing workflows with the following tutorials:

- [Constructing a simple worfklow](/tutorials/simple)
- [Running workflows](/tutorials/runningworkflows)
- [Building a bioinformatics workflow](/tutorials/alignsortedbam)

