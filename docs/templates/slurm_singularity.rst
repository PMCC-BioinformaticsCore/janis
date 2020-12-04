Slurm Singularity
=================

Template ID: ``slurm_singularity``



Quickstart
-----------

Take note below how to configure the template. This quickstart only includes the fields you absolutely require. This will write a new configuration to ``~/.janis.conf``. See `configuring janis <https://janis.readthedocs.io/en/latest/references/configuration.html>`__ for more information.

.. code-block:: bash

   janis init slurm_singularity
   
   # or to find out more information
   janis init slurm_singularity --help

OR you can insert the following lines into your template:

.. code-block:: yaml

   template:
     id: slurm_singularity



Fields
-------



**Optional**

=============================  ===================================  ==========================================  ========================================================================================================================================
ID                             Type                                 Default                                     Documentation
=============================  ===================================  ==========================================  ========================================================================================================================================
container_dir                  <class 'str'>                                                                    Location where to save and execute containers to, this will also look at the env variables 'CWL_SINGULARITY_CACHE', 'SINGULARITY_TMPDIR'
intermediate_execution_dir     <class 'str'>
queues                         typing.Union[str, typing.List[str]]                                              A single or list of queues that work should be submitted to
mail_program                                                                                                    Mail program to pipe email to, eg: 'sendmail -t'
send_job_emails                <class 'bool'>                       False                                       (requires JanisConfiguration.notifications.email to be set) Send emails for mail types END
catch_slurm_errors             <class 'bool'>                       True                                        Catch Slurm errors (like OOM or walltime)
build_instructions             <class 'str'>                        singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
singularity_load_instructions                                                                                   Ensure singularity with this command executed in shell
max_cores                                                                                                       Maximum number of cores a task can request
max_ram                                                                                                         Maximum amount of ram (GB) that a task can request
max_duration                                                                                                    Maximum amount of time in seconds (s) that a task can request
can_run_in_foreground          <class 'bool'>                       True
run_in_background              <class 'bool'>                       False
sbatch                         <class 'str'>                        sbatch                                      Override the sbatch command
=============================  ===================================  ==========================================  ========================================================================================================================================

