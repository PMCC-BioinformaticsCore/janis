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

==================  ===================================  ====================================================  ===============
ID                  Type                                 Default                                               Documentation
==================  ===================================  ====================================================  ===============
executionDir        <class 'str'>
queues              typing.Union[str, typing.List[str]]  prod_med,prod
containerDir        <class 'str'>                        /config/binaries/singularity/containers_devel/janis/
singularityVersion  <class 'bool'>                       3.4.0
catchSlurmErrors    <class 'bool'>                       True
sendSlurmEmails     <class 'bool'>                       False
max_workflow_time   <class 'int'>                        14400
==================  ===================================  ====================================================  ===============

