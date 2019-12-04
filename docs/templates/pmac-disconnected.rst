Pmac-Disconnected
=================

Template ID ``pmac-disconnected``

Fields
-------

**Required**

============  =============  ===============
ID            Type           Documentation
============  =============  ===============
executionDir  <class 'str'>
============  =============  ===============

**Optional**

==================  ===================================  ====================================================  ===============
ID                  Type                                 Default                                               Documentation
==================  ===================================  ====================================================  ===============
queues              typing.Union[str, typing.List[str]]  prod_med,prod
containerDir        <class 'str'>                        /config/binaries/singularity/containers_devel/janis/
singularityVersion  <class 'str'>                        3.4.0
catchSlurmErrors    <class 'bool'>                       True
sendSlurmEmails     <class 'bool'>                       False
==================  ===================================  ====================================================  ===============

