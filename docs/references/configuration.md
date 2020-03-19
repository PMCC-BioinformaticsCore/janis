# Configuring Janis

When we talk about configuring Janis, we're really talking about how to configure the [Janis assistant](https://github.com/PMCC-BioinformaticsCore/janis-assistant) and hence Cromwell / CWLTool to interact with batch systems (eg: Slurm / PBS Torque) and container environments (Docker / Singularity).

Janis has built in templates for the following compute environments:

- Local (Docker or Singularity)
    - The only environment compatible with CWLTool.
- Slurm (Singularity only)
- PBS / Torque (Singularity only)

In an extension of Janis ([janis-templates](https://github.com/PMCC-BioinformaticsCore/janis-templates)), we've produced a number of location specific templates, usually with sensible defaults.

See this list for a full list of templates: [https://janis.readthedocs.io/en/latest/templates/index.html](https://janis.readthedocs.io/en/latest/templates/index.html)

## Location of config


By default,
 
- Your configuration directory is placed at `~/.janis/`.
- Your configuration path is a file called `janis.conf` in this directory, eg: `~/.janis/janis.conf`  


Both `janis run / translate` allow you to provide a location to a config using the `-c / --config` parameter. 

In addition, you can also configure the following environment variables:


- `JANIS_CONFIGPATH` - simply the path to your config file
- `JANIS_CONFIGDIR` - The configuration path is determined by `$JANIS_CONFIGDIR/janis.conf`



## Initialising a template




## Other Janis options

This guide has a number of subsections for:

- `environment`
- `notifications`
- `cromwell`
- `recipes`
- `call_caching`
- `template` (previously discussed)


There are a number of high level keys you can use to configure Janis:

- `config_dir`: (default: `~/.janis/`) (equivalent env: `JANIS_CONFIGDIR`) - places internal database, downloaded workflows and other information in this directory.
- `output_dir`: Instead of specifying `-o` on every run, specifying this option will use an output directory with the following structure: `$output_dir/<workflow name>/yyyymmdd_hhMMss_<wid>` (where `wid` is the workflow ID).
- `execution_dir`: Sometimes it's useful to perform the execution in a different directory to where you'd like to place the outputs. One reason might be that Cromwell will copy (rather than hard link) your input files if they don't exist on the same drive) 
- `engine`: (OPTIONS: `cromwell` / `cwltool`) (overridable with `--engine ENGINE`) Which engine to use by default.
- `run_in_background`: (overridable with the `--background` and `--foreground` optiopns) Should workflows run in the background by default.


```yaml
config_dir: <default: ~/.janis/>
output_dir: <default output directory prefix>
execution_dir: <
engine: <cromwell / cwltool>

run_in_background: <boolean>
```


#### `output_dir`: Default output directory



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


### Existing Cromwell instance

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

## Call caching

## Recipes / Common inputs

Often between a number of different workflows, you have a set of inputs that you want to apply multiple times.
