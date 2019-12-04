Spartan-Disconnected
====================

Template ID ``spartan-disconnected``

Fields
-------

**Required**

============  ==============  ================================
ID            Type            Documentation
============  ==============  ================================
executionDir  <class 'type'>  execution directory for Cromwell
============  ==============  ================================

**Optional**

==================  ===================================  ====================================================  ========================================================
ID                  Type                                 Default                                               Documentation
==================  ===================================  ====================================================  ========================================================
queues              typing.Union[str, typing.List[str]]  physical                                              The queue to submit jobs to
containerDir        <class 'str'>                        /config/binaries/singularity/containers_devel/janis/
singularityVersion  <class 'str'>                        3.2.0-spartan_gcc-6.2.0
sendSlurmEmails     <class 'bool'>                       True                                                  Send SLURM job emails to the listed email address
catchSlurmErrors    <class 'bool'>                       True                                                  Fail the task if Slurm kills the job (eg: memory / time)
max_cores           <class 'int'>                        32                                                    Override maximum number of cores (default: 32)
max_ram             <class 'int'>                        508                                                   Override maximum ram (default 508 [GB])
==================  ===================================  ====================================================  ========================================================

