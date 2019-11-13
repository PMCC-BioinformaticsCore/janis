# Containerising a tool

Portability and reproducibility are important aspects of `janis`. Building and using containers for your tool is a great way to ensure compatibility and portability across institutes.

## Getting started

In this tutorial, we're going to introduce containers and go through different ways to find or build containers for your tool. 

### What you'll need

You'll need an installation of Docker, and a Dockerhub account for uploading your tool to the cloud. Other cloud-container vendors such as [quay.io](https://quay.io/) may be used instead of Dockerhub. 

You'll also need your tool, the dependencies and know how to install it.

### Licenses and warnings

You should have all relevant permission to be able to distribute the tool that you are wrapping. Most [open source licenses](https://choosealicense.com/appendix/#distribution) allow you distribute a compiled or derived version, but may have some restrictions.

### Versioning

It's important that we consider the versions of our tool when we construct our container. A large factor of this project is reproducibility, and this is only possible if we use the same tools. We highly discourage the use of the `latest` tag, and instead tag the container with its explicit version.

> Note: The same docker tag doesn't ensure you get the same container each time as a user can upload a new build for the same tag. Ideally you could use the Docker [`digest`](https://success.docker.com/article/images-tagging-vs-digests) however not all container software support this.

## Let's get started!

### What are containers
 
 Containers are allow an operating system and software to be virtualised on a different computer (the host). Although containers can be used for a variety of purposes, we use them as a "sandbox" for us to run analysis in. The important aspect for us is that we can guarantee reproducibility with the same container across compute platforms. 

In this guide we'll be creating Docker containers, however other virtualisation software (such as Singularity, uDocker or Shifter) can read this common ([OCI compliant](https://www.opencontainers.org/)) standard.

### How do I get a container?

When we look for a container, we want to make sure it's high quality and is as described. Some services (such as the community supported [Biocontainers](https://biocontainers.pro/#/registry)) produce high quality containers. 

Here is a number of criteria we're looking for in a container:

- From a reputable source.
- The Dockerfile is available (gives us trust that there's not malicious software in there).
- The versions available are sensible and match the software versions.
	- A container with only the `latest` tag does not satisfy this criteria.

Here is our recommendation for finding a container given the listed criteria:

1. Look for a container from the tool provider (maybe through their GitHub).
2. Look for a container from a reputable source on [Dockerhub](https://hub.docker.com/search/?q=&type=image) or [quay.io](https://quay.io/search?q=).
3. Find a third-party container that meets all the criteria above.

If you are unable to find a container, you might need to build one (and maybe find a tech savvy colleague to help you).


## Building a container

We're going to build a Docker container as at the moment they are the most portable option. Each Docker container has a recipe, called a [`Dockerfile`](https://docs.docker.com/engine/reference/builder/), it's a special file that tells Docker how a container is built, and how it should be executed. 

Each container should have a base, ie an operating system or even another docker container.

Here are the following types of containers we'll create:

- [Python module](#python-module)
- [Python script](#python-script)
- [Downloadable binary w/ requirements](#runnable-jar)
- [Makeable program](#makeable-program)

### Basic setup

For all of the following guides, you'll need to create a folder for your Dockerfile and other dependencies. It's important you place your dependencies within this directory (or subdirectories) to avoid bloating your [build context](https://medium.com/lucjuggery/docker-tips-about-the-build-context-dbc76505e178).

```bash
tooldir=mytool-docker
mkdir $tooldir && cd $tooldir
touch Dockerfile
```

### Python module

We're going to package up a Python module that you can run from the console. For this example we'll wrap Janis in a container.

> Within the `setup.py`, it has a dictionary value within the `entry_points` kwarg and a value for the `console_script` key. For example, Janis-assitant [setup.py:24](https://github.com/PMCC-BioinformaticsCore/janis-assistant/blob/master/setup.py#L24 ) looks like this:
> ```python
> entry_points={
>     "console_scripts": ["janis=janis_assistant.cli:process_args"],
> }
> ```

Given our application is compatible with Python 3.7, we simply `pip install` our module on top, place this in our `Dockerfile`:
```docker
FROM python:3.7
RUN pip install janis-pipelines
CMD janis
```

> Although `janis-pipelines` is compatible with the alpine Python 3.7 container, we'll use the full version, see "[The best Docker base image for your Python application](https://pythonspeed.com/articles/base-image-python-docker-images/)" for more information.

We can build this container with the following command which will
- Use the build context: `.` 
- _Automatically_ build the file called `Dockerfile` relative to the build context
- Give the container the tag: `yourname/packagename` 

```bash
docker build -t yourname/janis-pipelines .
```

We can test the container works by running:
```bash
docker run yourname/janis-pipelines janis -v
#--------------------  ------
#janis-core            v0.7.1
#janis-assistant       v0.7.7
#janis-unix            v0.7.0
#janis-bioinformatics  v0.7.1
#--------------------  ------
```

OR:
```bash
docker run -it --entrypoint /bin/sh yourname/janis-pipelines
## inside the container
$ janis translate hello wdl 
# translation here
$ exit # exit the container
```



### Python script

We have a single python script that we want to run inside a Docker container. This time we're going to add our files to the `/opt` directory (in the container), add that directory to the `$PATH` variable and then we can simply call your docker container.

We're going to create a little Python program to write something to the console and then end. Save the following Python program (including the [Shebang](https://en.wikipedia.org/wiki/Shebang_(Unix))) into the directory with our Dockerfile.

**`hello.py`**:
```python
#!/usr/bin/env python3
print("Hello, World")
```

We need to make `hello.py` executable by:
```bash
chmod +x hello.py
```
> We can test this worked by calling our script, by running `./hello.py` in our terminal

We'll now edit our `Dockerfile` with the following steps:
- Use a `python:3.7` base.
- Add the file into container.
- Modify the path to include the script directory (`/opt`).
- Use the `hello.py` script an entry-point.

```docker
FROM python:3.7.5-alpine
ADD hello.py /opt/
ENV PATH="/opt:${PATH}"
CMD hello.py
```

We can build this container with the following command:

> Don't forget the full stop `.` at the end of the docker build command

```bash
docker build -t yourname/hello .
```

You can test that it worked in two ways:

1. Run the container (using the entrypoint (CMD) `hello.py`):
```bash
docker run yourname/hello
# Hello, World!
```

2. Go into the container, and run the script:
```bash
docker run -it --entrypoint /bin/sh yourname/hello
## inside the container
$ hello.py
# Hello, World!
$ exit # exit the container
```

### Downloadable Binary with Requirements

> This section is still under construction.

This style of container is very similar to the Python script. We'll `ADD` our binary (potentially from the web) to `/opt/toolname`, add this to the path and add an entry point.

You will need to consider which Base OS to use, for now we'll use alpine as it's very lightweight, but you might need to consider `ubuntu` or `centos`. We'll assume that your tool `mytool` runs out of a `bin` folder.

**`Dockerfile`**

```docker
FROM alpine:latest
RUN mkdir -p /opt/mytool/
ADD https://mytool.net/releases/tool.zip /opt/
ENV PATH="/opt/mytool/bin:${PATH}"	# assume ./mytool/bin/
CMD mytool
```

### Makeable Program

> This section is under construction, please refer to an example [Samtools Docker](https://github.com/PMCC-BioinformaticsCore/scripts/blob/master/dockers/samtools/Dockerfile) for more information:

**`Dockerfile`**

```docker
FROM ubuntu:16.04

MAINTAINER Michael Franklin <michael.franklin@petermac.org>

RUN apt-get update -qq \
  && apt-get install -qq bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev

ENV SAMTOOLS_VERSION 1.9

LABEL \
  version="${SAMTOOLS_VERSION}" \
  description="Samtools image for use in Workflows"
  
RUN cd /opt/ \
	&& wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
	&& cd samtools-${SAMTOOLS_VERSION}/ \
	&& make && make install

ENV PATH="/opt/samtools-${SAMTOOLS_VERSION}/:${PATH}"
```


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


### The less steps your run, the smaller the container

By joining `RUN` commands together using the `&&` operator, your total containers may be smaller in size.

