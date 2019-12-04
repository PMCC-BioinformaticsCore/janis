Pmac
====

Template ID ``pmac``

Fields
-------

**Required**

============  =============  ===================
ID            Type           Documentation
============  =============  ===================
executionDir  <class 'str'>  Execution directory
============  =============  ===================

**Optional**

============================  ===================================  ====================================================  ======================================================================
ID                            Type                                 Default                                               Documentation
============================  ===================================  ====================================================  ======================================================================
queues                        typing.Union[str, typing.List[str]]  prod_med,prod                                         The queue to submit jobs to
containerDir                  <class 'str'>                        /config/binaries/singularity/containers_devel/janis/  [OPTIONAL] Override the directory singularity containers are stored in
singularityVersion            <class 'str'>                        3.4.0                                                 The version of Singularity to use on the cluster
sendSlurmEmails               <class 'bool'>                       False                                                 Send Slurm job notifications using the provided email
catchSlurmErrors              <class 'bool'>                       True                                                  Fail the task if Slurm kills the job (eg: memory / time)
singularityBuildInstructions                                                                                             Sensible default for PeterMac template
max_cores                     <class 'int'>                        40                                                    Override maximum number of cores (default: 32)
max_ram                       <class 'int'>                        256                                                   Override maximum ram (default 508 [GB])
============================  ===================================  ====================================================  ======================================================================

