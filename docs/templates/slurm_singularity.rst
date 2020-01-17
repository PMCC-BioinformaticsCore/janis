Slurm_Singularity
=================

Template ID ``slurm_singularity``

Fields
-------

**Required**

=============  =============  ==================================================
ID             Type           Documentation
=============  =============  ==================================================
container_dir  <class 'str'>  Location where to save and execute containers from
=============  =============  ==================================================

**Optional**

=============================  ===================================  ==========================================  ==========================================================================================
ID                             Type                                 Default                                     Documentation
=============================  ===================================  ==========================================  ==========================================================================================
execution_dir                  <class 'str'>
queues                         typing.Union[str, typing.List[str]]                                              A single or list of queues that work should be submitted to
mail_program                                                                                                    Mail program to pipe email to, eg: 'sendmail -t'
send_job_emails                <class 'bool'>                       True                                        (requires JanisConfiguration.notifications.email to be set) Send emails for mail types END
catch_slurm_errors             <class 'bool'>                       False                                       Catch Slurm errors (like OOM or walltime)
build_instructions             <class 'str'>                        singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
singularity_load_instructions                                                                                   Ensure singularity with this command executed in shell
limit_resources                <class 'bool'>                       False                                       Limit resources with singularity using cgroups (REQUIRES ROOT)
max_cores                                                                                                       Maximum number of cores a task can request
max_ram                                                                                                         Maximum amount of ram (GB) that a task can request
=============================  ===================================  ==========================================  ==========================================================================================

