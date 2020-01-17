Singularity
===========

Template ID ``singularity``

Fields
-------

**Required**

=============  ==============  ===============
ID             Type            Documentation
=============  ==============  ===============
container_dir  <class 'type'>
=============  ==============  ===============

**Optional**

===========================  =============  ==========================================  ==================================================================================
ID                           Type           Default                                     Documentation
===========================  =============  ==========================================  ==================================================================================
singularityLoadInstructions                                                             Ensure singularity with this command executed in shell
containerBuildInstructions   <class 'str'>  singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
mailProgram                  <class 'str'>                                              Mail program to pipe email to, eg: 'sendmail -t'
===========================  =============  ==========================================  ==================================================================================

