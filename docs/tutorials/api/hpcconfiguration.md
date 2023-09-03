# Configuring Janis for HPCs

The Janis assistant currently implements Cromwell as an execution engine which supports HPCs.

The best way to use Janis and your HPC is to use one of the provided configurations. Currently, all of these configs use `singularity` to manage docker containers.

See a list of templates [here](https://janis.readthedocs.io/en/latest/templates/index.html).

## Example: Slurm

- Template ID: `slurm_singularity`
- Documentation: `https://janis.readthedocs.io/en/latest/templates/slurm_singularity.html`
- CLI help: `janis init slurm_singularity --help`.

Example to configure this:
```bash
# This will write a janis.conf to your $HOME/.janis/janis.conf
$ janis init slurm_singularity

$ cat ~/.janis/janis.conf
# engine: cromwell
# notifications:
#   email: null
# template:
#   catch_slurm_errors: true
#   id: slurm_singularity
#   max_workflow_time: 20100
#   sbatch: sbatch
#   send_job_emails: false
```


### Background mode

The slurm_singularity template constructs 