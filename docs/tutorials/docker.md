# Dockerizing a tool

Portability and reproducibility are very important aspects of `janis`, dockerizing your tool is a great way to ensure compatibility and remove the burden of modules.

## Getting started

In this tutorial, we're going to discuss Docker and work through the process of creating a docker container of some different types of tools.

### What you'll need

You'll need an installation of Docker, and a Dockerhub account for uploading your tool to the cloud. Other cloud-container vendors such as [quay.io](https://quay.io/) may be used instead of Dockerhub. It helps to have some basic knowledge

You'll also need, a tool, its dependency requirements and know how to install it.

### Licenses and warnings

You should have all relevant permission to be able to distribute the tool that you are wrapping. Most [open source licenses](https://choosealicense.com/appendix/#distribution) allow you distribute a compiled or derived version, but may have some restrictions.

### Versioning

It's important that we're particular with versioning of our tool when we construct our container. A large factor of this project is reproducibility, and this is only possible if we use the same tools. We highly recommend that you avoid the `latest` tag, and instead explicitly tag the container.

> Note: The same docker tag doesn't ensure you get the same container each time as a user can upload a new build for the same tag. Ideally you could use the Docker [`digest`](https://success.docker.com/article/images-tagging-vs-digests) however not all programs support this.

## Let's get started!

Each Docker container has a recipe, called a [`Dockerfile`](https://docs.docker.com/engine/reference/builder/), it's a special file that tells Docker how a container is built, and how it should be executed. 

Each container should have a base, ie an operating system or even another docker container.

Here are the following types of Docker containers we'll create:
- [Python module](#python-module)
- [Python script](#python-script)
- [Downloadable binary w/ requirements](#runnable-jar)
- [Makeable program](#makeable-program)

### Python module

We're going to package up a Python module that you can run from the console, hence within the `setup.py`, it has a dictionary value within the `entry_points` kwarg and a value for the `console_script` key. For example, CWLTool's [setup.py:76](https://github.com/common-workflow-language/cwltool/blob/master/setup.py#L76) looks like this:
```python
entry_points={
    'console_scripts': ["cwltool=cwltool.main:run"]
},
```

We need to choose a base for our Docker container, and so given our application is compatible with Python 3.7, then we simply `pip install` our module on top. The following Dockerfile should be saved with the name: `Dockerfile`:
```docker
FROM python:3.7
RUN pip install package-name
```

We can then build this container with the following command which will
- Use the build context: `.` 
- _Automatically_ build the file called `Dockerfile` relative to the build context
- Give the container the tag: `yourname/packagename` 

```bash
docker build -t yourname/packagename .
```

### Python script

We have a single python script that we want to run inside a Docker container. This time we're going to copy files into a specific directory, add that directory to the `$PATH` variable and then we can simply call your docker container.

We're going to create a little Python program to write something to the console and then end. Save the following Python program (including the [Shebang](https://en.wikipedia.org/wiki/Shebang_(Unix))) into an empty directory as `hello.py`:

```python
#!/usr/bin/env python3
print("Hello, World")
```

We need to make `hello.py` executable by:
```bash
chmod +x hello.py
```
> We can test this worked by calling our script, by running `./hello.py` in our terminal

We'll now create a `Dockerfile` with a `python:3.7` base, copy the file into container, add the script to path and run it as an entry-point.

```docker
FROM python:3.7
RUN mkdir /install_dir
ENV PATH="install_dir:${PATH}"

COPY hello.py /install_dir/hello.py
CMD hello.py
```

We can build this container with the following command:
```bash
docker build -t yourname/hello .
```

### Downloadable Binary with Requirements

_This section is under construction_


### Makeable Program

_This section is under construction_


## Additional helpful hints

### Keeping Dockerfiles separate from source code

You can keep your Dockerfile in a different place to your file (with restrictions) however you must change a few arguments:

- Your build context (the last argument) should be in a place that can access both the Dockerfile and program files as subdirectories
- The `COPY hello.py` command in the Dockerfile must be changed to be relative to the build context.
- You can then supply `-f` argument to point to the Dockerfile, relative to the build context:

Example: If the python `hello.py` file is stored in a subdirectory called `src`, and the Dockerfile is stored in a subdirectory called `docker_stuff`, you could use the following:

Dockerfile:
> ```
> COPY src/hello.py /install_dir/hello.py
> ```

Build:
> Change into the parent directory of `src/` and `docker_stuff/`
> ```bash
> docker build -t yourname/hello -f dockerstuff/Dockerfile .
> ```
