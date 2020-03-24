# Configuring Janis Assistant

When we talk about configuring Janis, we're really talking about how to configure the [Janis assistant](https://github.com/PMCC-BioinformaticsCore/janis-assistant) and hence Cromwell / CWLTool to interact with batch systems (eg: Slurm / PBS Torque) and container environments (Docker / Singularity).

Janis has built in templates for the following compute environments:

- Local (Docker or Singularity)
    - The only environment compatible with CWLTool.
- Slurm (Singularity only)
- PBS / Torque (Singularity only)

In an extension of Janis ([janis-templates](https://github.com/PMCC-BioinformaticsCore/janis-templates)), we've produced a number of location specific templates, usually with sensible defaults.

See this list for a full list of templates: [https://janis.readthedocs.io/en/latest/templates/index.html](https://janis.readthedocs.io/en/latest/templates/index.html)

## Syntax

The config should be in YAML, and can contain nested dictionaries.

**NOTE**: In this guide we might use dots (`.`) to refer to a nested key. For example, `notifications.email` refers to the following structure:

```yaml
notifications:
  email: <value>
```

## Location of config


By default,
 
- Your configuration directory is placed at `~/.janis/`.
- Your configuration path is a file called `janis.conf` in this directory, eg: `~/.janis/janis.conf`  


Both `janis run / translate` allow you to provide a location to a config using the `-c / --config` parameter. 

In addition, you can also configure the following environment variables:


- `JANIS_CONFIGPATH` - simply the path to your config file
- `JANIS_CONFIGDIR` - The configuration path is determined by `$JANIS_CONFIGDIR/janis.conf`



## Initialising a template

```bash
janis init <template> [...options]
```

By default this will place our config at `~/.janis/janis.conf`. If there's already a config there, it will NOT override it.
 
 - `-o` / `--output`: output path to write to, default (`~/.janis/janis.conf`).
 - `-f` / `--force`: Overwrite the config if one exists at the output path.
 - `--stdout`: Write the output to stdout. 

### Example Slurm + Singularity

Let's choose the [`slurm_singularity`](https://janis.readthedocs.io/en/latest/templates/slurm_singularity.html) template. We can look at it's signature via the command line:

```bash
janis init slurm_singularity -h

# required arguments:
#   --container_dir CONTAINER_DIR
#                         Location where to save and execute containers from
# 
# optional arguments:
#   --execution_dir EXECUTION_DIR
#   --queues QUEUES       A single or list of queues that work should be
#                         submitted to
#   --mail_program MAIL_PROGRAM
#                         Mail program to pipe email to, eg: 'sendmail -t'
#   --send_job_emails     (default: True) (requires
#                         JanisConfiguration.notifications.email to be set) Send
#                         emails for mail types END
#   --catch_slurm_errors  (default: False) Catch Slurm errors (like OOM or
#                         walltime)
#   --build_instructions BUILD_INSTRUCTIONS
#                         (default: singularity pull $image docker://${docker})
#                         Instructions for building singularity, it's
#                         recommended to not touch this setting.
#   --singularity_load_instructions SINGULARITY_LOAD_INSTRUCTIONS
#                         Ensure singularity with this command executed in shell
#   --max_cores MAX_CORES
#                         Maximum number of cores a task can request
#   --max_ram MAX_RAM     Maximum amount of ram (GB) that a task can request
#   --can_run_in_foreground
#                         (default: False) None
#   --run_in_background   (default: True) None
```

We _MUST_ pass a value to `--container_dir` (which is where our singularity containers will be downloaded to), with the other parameters being optional.

We could then initialise our template with:

```bash
janis init slurm_singularity --container_dir /shared/path/to/singularity_containers/
```


## Other Janis options

This guide has a number of subsections for:

- `environment`
- `notifications`
- `cromwell`
- `call_caching`
- `recipes`
- `template` (previously discussed)


There are a number of high level keys you can use to configure Janis:

- `config_dir`: (default: `~/.janis/`) (equivalent env: `JANIS_CONFIGDIR`) - places internal database, downloaded workflows and other information in this directory.
- `output_dir`: Instead of specifying `-o` on every run, specifying this option will use an output directory with the following structure: `$output_dir/<workflow name>/yyyymmdd_hhMMss_<wid>` (where `wid` is the workflow ID).
- `execution_dir`: By default, your execution is placed in `<output dir>/janis/execution`, however sometimes it's useful to perform the execution in a different directory. One reason might be that Cromwell will copy (rather than hard link) your input files if they don't exist on the same drive) 
- `engine`: (OPTIONS: `cromwell` / `cwltool`) (overridable with `--engine ENGINE`) Which engine to use by default.
- `run_in_background`: (overridable with the `--background` and `--foreground` options) Ask janis to submit workflows in the background.

For example:

```yaml
config_dir: <default: ~/.janis/>
output_dir: <default output directory prefix>
execution_dir: <execution dir override>
engine: <cromwell / cwltool>
run_in_background: <boolean>
```

### Environment variables

> Non exhaustive list, some additional options are in the subsections below.

- `JANIS_EXCECUTIONDIR`: equivalent to `execution_dir` as above.
- `JANIS_OUTPUTDIR`: equivalent to `output_dir` as above.

```bash
export JANIS_EXCECUTIONDIR=/path/to/alternative/execution/
export JANIS_OUTPUTDIR=/path/to/default/output/
```


## Environment

> These options will sit under a `environment` key, see the example for more information.

There are a couple of environment options you might want to alter:

- `max_cores`: Limit the maximum amount of requested cores. 
- max_ram

```yaml
# other params
environment:
  max_cores: int
  max_ram: int

```



## Notifications

- `email`: An email address to send notifications to
- `mail_program`: Janis constructs the entire email (including the headers), and passes it through `stdin` to this program. Most cluster have `sendmail -t` installed, the `-t` parameter is important to "Extract recipients from message headers". 


Example:

```yaml
notifications: 
  email: emailtosendupdates@domain.com
  mail_program: sendmail -t
```

## Cromwell

Sometimes Cromwell can be hard to configure from a simple template, so we've exposed some extra common options. Feel free to [raise an issue](https://github.com/PMCC-BioinformaticsCore/janis-assistant/issues/new) if you have questions or ideas.

- `jar`: A location to a specific Cromwell Jar. By default, it looks for jars in `JANIS_CONFIGDIR` (default: ~/.janis/`), and finds the highest version. 
- `config_path`: Path to an existing Cromwell Configuration. You might have a configuration file you've written (and that the templates here don't cover very well). Potentially a preconfigured version for the cloud. 
- `url`: (URL + PORT) Or maybe you've already got an existing Cromwell instance (running on the cloud as well) but you'd still like the benefits that Janis provides (auto generation, submission, progress tracking and output copying). 
- `memory_mb`: Increase the default heap size (`-xmx`) of Cromwell. The minimum (`-xms`) is raised to half of this value. (See it's usage here: [engines/cromwell/main.py#L163-L166](https://github.com/PMCC-BioinformaticsCore/janis-assistant/blob/edb2044705f8221f52e2520eafc0779efcbf286a/janis_assistant/engines/cromwell/main.py#L163-L166))
- `call_caching_method`: (Default: `file`, but will become [`fingerprint`](https://github.com/broadinstitute/cromwell/pull/5450)) ([Cromwell Docs](https://cromwell.readthedocs.io/en/stable/Configuring/#local-filesystem-options)).


```yaml
cromwell:
  jar: <path-to-cromwell.jar>
  config_path: <path-to-cromwell.conf>
  url: localhost:8000
  memory_mb: 1500 # don't include MB or anything with it, this MUST be an int
  call_caching_method: 
```


### Existing Cromwell instance

In addition to the previous `cromwell.url` method, you can also manage this through the command line.
To configure Janis to submit to an existing Cromwell instance (eg: one already set up for the cloud), the CLI has a mechanism for setting the `cromwell_url`:

```bash
urlwithport="127.0.0.1:8000"
janis run --engine cromwell --cromwell-url $urlwithport hello
```

OR

**`~/janis.conf`**

```yaml
engine: cromwell
cromwell:
  url: 127.0.0.1:8000
```

### Overriding Cromwell JAR

In additional to the previous `cromwell.jar`, you can set the location of the Cromwell JAR through the environment variable `JANIS_CROMWELLJAR`:

```bash
export JANIS_CROMWELLJAR=/path/to/cromwell.jar
```

## Call caching

Please refer to this guide for configuring Janis + (Cromwell / CWLTool) to use call-caching: [Call Caching](https://janis.readthedocs.io/en/latest/references/callcaching.html).

## Recipes / Common inputs

Often between a number of different workflows, you have a set of inputs that you want to apply multiple times. For that, we have `recipes`. 

> Recipes only match on the input name. Your input names _MUST_ be consistent through your pipelines for this concept to be useful.

For example, everytime I run the [`WGSGermlineGATK`](https://janis.readthedocs.io/en/latest/pipelines/wgsgermlinegatk.html) pipeline with `hg38`, I know I want to provide the same reference files.

There are a few ways to configure this:

### Directories

You can setup a directory with a series of `*.yaml` files, tell Janis about it; Janis will look up the file names in that directory.

For the following directory: 

```
recipe_directory/
    hg38.yaml
    hg19.yaml
    mm10.yaml
```

You could then say `janis run -r hg38 [...options] WGSGermlineGATK`.

Environment variable:

```bash
export JANIS_RECIPEDIRECTORY=/path/to/recipe_directory/
```

OR, in your `janis.conf`

```yaml
recipes:
  directories:
  - /path/to/recipe_directory/
```

### Paths

Create a file called `recipes.yaml` with the following structure:

```yaml
hg38:
  reference: /path/to/hg38/reference.fasta
  type: hg38
hg19:
  reference: /path/to/hg19/reference.fasta
  type: hg19
```

Environment variable:

```bash
export JANIS_RECIPEPATHS="comma,separated,/path/to/recipes.yaml
```

OR, in your `janis.conf`

```yaml
recipes:
  paths:
  - /path/to/recipes.yaml 
```

### Directly in your `janis.conf`

```yaml
recipes:
  recipes:
    hg38:
      reference: /path/to/hg38/reference.fasta
      type: hg38
    hg19:
      reference: /path/to/hg19/reference.fasta
      type: hg19
```






