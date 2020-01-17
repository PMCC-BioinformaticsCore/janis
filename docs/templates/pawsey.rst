Pawsey
======

Template ID ``pawsey``

Fields
-------

**Required**

=============  =============  ==================================================
ID             Type           Documentation
=============  =============  ==================================================
container_dir  <class 'str'>  Location where to save and execute containers from
=============  =============  ==================================================

**Optional**

==============================  ===================================  ==========================================  ==========================================================================================
ID                              Type                                 Default                                     Documentation
==============================  ===================================  ==========================================  ==========================================================================================
execution_dir                   <class 'str'>
queues                          typing.Union[str, typing.List[str]]  workq                                       A single or list of queues that woork should be submitted to
singularity_version             <class 'str'>                        3.3.0                                       Version of singularity to load
catch_slurm_errors              <class 'bool'>                       True                                        Catch Slurm errors (like OOM or walltime)
send_job_emails                 <class 'bool'>                       True                                        (requires JanisConfiguration.notifications.email to be set) Send emails for mail types END
singularity_build_instructions  <class 'str'>                        singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
max_cores                       <class 'int'>                        28                                          Maximum number of cores a task can request
max_ram                         <class 'int'>                        128                                         Maximum amount of ram (GB) that a task can request
==============================  ===================================  ==========================================  ==========================================================================================

