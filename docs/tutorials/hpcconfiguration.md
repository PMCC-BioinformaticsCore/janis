# Configuring Janis for HPCs

The Janis assistant currently implements Cromwell as an execution engine which supports HPCs.

The best way to use Janis and your HPC is to use one of the provided configurations. Currently, all of these configs use `singularity` to manage docker containers.

See a list of templates [here](https://janis.readthedocs.io/en/latest/templates/index.html).

## Example: Slurm

- Template ID: `slurm_singularity`
- Documentation: `https://janis.readthedocs.io/en/latest/templates/slurm_singularity.html`
- CLI help: `janis init slurm_singularity --help`.

Minimum configuration:

- `executionDir`: Where should the execution take place
- `containerDir`: Where should singularity containers be stored.

Example to configure this:
```bash
janis init slurm_singularity --executionDir /path/to/executionDir --containerDir /path/to/containerDir
```

More configuration:

- `queues`: A queue or list of queues to submit to
- `max_cores`: The maximum amount of acores 
- `mail_program`: Mail program for which an email is piped to (eg: `sendmail -t`)
- `buildInstructions`: How to build a singularity container, default: `singularity pull $image docker://${docker}`. Don't change this unless you absolutely have it.

