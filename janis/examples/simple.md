# Simple Workflow

In this tutorial

### Designing the workflow

A simple untar workflow will give us an understanding of how inputs, steps, outputs, types and connections work.
The order of the workflow should be:
1. [Untar](https://janis.readthedocs.io/en/latest/tools/unix/untar.html) the file.
2. [Compile](https://janis.readthedocs.io/en/latest/tools/unix/javacompiler.html) the resulting files.
3. [Tar](https://janis.readthedocs.io/en/latest/tools/unix/tar.html) the compiled files.

I've directly linked the tools to their documentation pages. 


### Using the documentation

The documentation is the best way to determine how the tools work. When tools are _wrapped_ (written), 
the author can declare a description, links to more information and more, the documentation will be generated
from these tool wrappers.

The documentation 


### Declaring the workflow

The first step is to import Janis:
```python
import janis as j
```



We declare the inputs

The inputs, steps and outputs are all declared in the same way.