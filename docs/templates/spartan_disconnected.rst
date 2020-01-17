Spartan_Disconnected
====================

Template ID ``spartan_disconnected``

Fields
-------

**Required**

=============  =============  ===============
ID             Type           Documentation
=============  =============  ===============
container_dir  <class 'str'>
=============  =============  ===============

**Optional**

===================  ===================================  =======================  ========================================================
ID                   Type                                 Default                  Documentation
===================  ===================================  =======================  ========================================================
execution_dir        <class 'str'>                                                 execution directory for Cromwell
queues               typing.Union[str, typing.List[str]]  physical                 The queue to submit jobs to
singularity_version  <class 'str'>                        3.2.0-spartan_gcc-6.2.0
send_job_emails      <class 'bool'>                       True                     Send SLURM job emails to the listed email address
catch_slurm_errors   <class 'bool'>                       True                     Fail the task if Slurm kills the job (eg: memory / time)
max_cores            <class 'int'>                        32                       Override maximum number of cores (default: 32)
max_ram              <class 'int'>                        508                      Override maximum ram (default 508 [GB])
===================  ===================================  =======================  ========================================================

