Pmac_Disconnected
=================

Template ID ``pmac_disconnected``

Fields
-------

**Required**

====  ======  ===============
ID    Type    Documentation
====  ======  ===============
====  ======  ===============

**Optional**

===================  ===================================  ====================================================  ===============
ID                   Type                                 Default                                               Documentation
===================  ===================================  ====================================================  ===============
execution_dir        <class 'str'>
queues               typing.Union[str, typing.List[str]]  prod_med,prod
container_dir        <class 'str'>                        /config/binaries/singularity/containers_devel/janis/
singularity_version  <class 'bool'>                       3.4.0
catch_slurm_errors   <class 'bool'>                       True
send_job_emails      <class 'bool'>                       False
max_cores            <class 'int'>                        40
max_ram              <class 'int'>                        256
max_workflow_time    <class 'int'>                        14400
===================  ===================================  ====================================================  ===============

