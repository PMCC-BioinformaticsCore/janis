Singularity
===========

Template ID ``singularity``

Fields
-------

**Required**

====  ======  ===============
ID    Type    Documentation
====  ======  ===============
====  ======  ===============

**Optional**

===========================  =============  ====================================================  ==================================================================================
ID                           Type           Default                                               Documentation
===========================  =============  ====================================================  ==================================================================================
executionDir
containerDir                 <class 'str'>  /config/binaries/singularity/containers_devel/janis/  Location where to save and execute containers from
singularityLoadInstructions                                                                       Ensure singularity with this command executed in shell
containerBuildInstructions   <class 'str'>  singularity pull $image docker://${{docker}}          Instructions for building singularity, it's recommended to not touch this setting.
mailProgram                  <class 'str'>                                                        Mail program to pipe email to, eg: 'sendmail -t'
===========================  =============  ====================================================  ==================================================================================

