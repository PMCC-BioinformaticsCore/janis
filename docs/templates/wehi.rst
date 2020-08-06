Wehi
====

Template ID: ``wehi``



Quickstart
-----------

Take note below how to configure the template. This quickstart only includes the fields you absolutely require. This will write a new configuration to ``~/.janis.conf``. See `configuring janis <https://janis.readthedocs.io/en/latest/references/configuration.html>`__ for more information.

.. code-block:: bash

   janis init wehi \
       --container_dir <value>
   
   # or to find out more information
   janis init wehi --help

OR you can insert the following lines into your template:

.. code-block:: yaml

   template:
     id: wehi
     container_dir: <value>



Fields
-------

**Required**

=============  =============  ===============
ID             Type           Documentation
=============  =============  ===============
container_dir  <class 'str'>
=============  =============  ===============

**Optional**

==========================  ===================================  =========  ==========================================================================================
ID                          Type                                 Default    Documentation
==========================  ===================================  =========  ==========================================================================================
intermediate_execution_dir  <class 'str'>                                   A location where the execution should take place
queues                      typing.Union[typing.List[str], str]             A single or list of queues that woork should be submitted to
singularity_version         <class 'str'>                        3.4.1      Version of singularity to load
catch_pbs_errors            <class 'bool'>                       True       Catch PBS errors (like OOM or walltime)
send_job_emails             <class 'bool'>                       True       (requires JanisConfiguration.notifications.email to be set) Send emails for mail types END
build_instructions                                                          Instructions for building singularity, it's recommended to not touch this setting.
max_cores                   <class 'int'>                        40         Maximum number of cores a task can request
max_ram                     <class 'int'>                        256        Maximum amount of ram (GB) that a task can request
==========================  ===================================  =========  ==========================================================================================

