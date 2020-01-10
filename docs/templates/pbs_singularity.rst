Pbs_Singularity
===============

Template ID ``pbs_singularity``

Fields
-------

**Required**

============  =============  ==================================================
ID            Type           Documentation
============  =============  ==================================================
containerDir  <class 'str'>  Location where to save and execute containers from
============  =============  ==================================================

**Optional**

===========================  ==============  ==========================================  ==========================================================================================
ID                           Type            Default                                     Documentation
===========================  ==============  ==========================================  ==========================================================================================
executionDir                 <class 'str'>
queue                        <class 'str'>                                               A queue that work should be submitted to
mail_program                                                                             Mail program to pipe email to, eg: 'sendmail -t'
sendJobEmails                <class 'bool'>  True                                        (requires JanisConfiguration.notifications.email to be set) Send emails for mail types END
catchPbsErrors               <class 'bool'>  True
buildInstructions            <class 'str'>   singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
singularityLoadInstructions                                                              Ensure singularity with this command executed in shell
max_cores                                                                                Maximum number of cores a task can request
max_ram                                                                                  Maximum amount of ram (GB) that a task can request
===========================  ==============  ==========================================  ==========================================================================================

