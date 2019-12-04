Pawsey-Disconnected
===================

Template ID ``pawsey-disconnected``

Fields
-------

**Required**

============  =============  ==================================================
ID            Type           Documentation
============  =============  ==================================================
executionDir  <class 'str'>
containerDir  <class 'str'>  Location where to save and execute containers from
============  =============  ==================================================

**Optional**

============================  ===================================  ==========================================  ==========================================================================================
ID                            Type                                 Default                                     Documentation
============================  ===================================  ==========================================  ==========================================================================================
queues                        typing.Union[str, typing.List[str]]  workq                                       A single or list of queues that woork should be submitted to
submissionQueue               <class 'str'>                        longq
singularityVersion            <class 'str'>                        3.3.0                                       Version of singularity to load
catchSlurmErrors              <class 'bool'>                       True                                        Catch Slurm errors (like OOM or walltime)
sendSlurmEmails               <class 'bool'>                       True                                        (requires JanisConfiguration.notifications.email to be set) Send emails for mail types END
singularityBuildInstructions  <class 'str'>                        singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
max_cores                     <class 'int'>                        28                                          Maximum number of cores a task can request
max_ram                       <class 'int'>                        128                                         Maximum amount of ram (GB) that a task can request
============================  ===================================  ==========================================  ==========================================================================================

