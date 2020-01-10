Spartan
=======

Template ID ``spartan``

Fields
-------

**Required**

============  =============  ===============
ID            Type           Documentation
============  =============  ===============
containerDir  <class 'str'>
============  =============  ===============

**Optional**

==================  ===================================  =======================  ========================================================
ID                  Type                                 Default                  Documentation
==================  ===================================  =======================  ========================================================
executionDir        <class 'str'>                                                 execution directory for Cromwell
queues              typing.Union[str, typing.List[str]]  physical                 The queue to submit jobs to
singularityVersion  <class 'str'>                        3.2.0-spartan_gcc-6.2.0
sendSlurmEmails     <class 'bool'>                       True                     Send SLURM job emails to the listed email address
catchSlurmErrors    <class 'bool'>                       True                     Fail the task if Slurm kills the job (eg: memory / time)
max_cores           <class 'int'>                        32                       Override maximum number of cores (default: 32)
max_ram             <class 'int'>                        508                      Override maximum ram (default 508 [GB])
==================  ===================================  =======================  ========================================================

