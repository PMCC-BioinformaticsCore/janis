Pawsey
======

Template ID: ``pawsey``


https://support.pawsey.org.au/documentation/display/US/Queue+Policies+and+Limits

Template for Pawsey. This submits Janis to the longq cluster. There is currently NO support
for workflows that run for longer than 4 days, though workflows can be resubmitted after this
job dies.

It's proposed that Janis assistant could resubmit itself


Quickstart
-----------

Take note below how to configure the template. This quickstart only includes the fields you absolutely require. This will write a new configuration to ``~/.janis.conf``. See `configuring janis <https://janis.readthedocs.io/en/latest/references/configuration.html>`__ for more information.

.. code-block:: bash

   janis init pawsey \
       --container_dir <value>
   
   # or to find out more information
   janis init pawsey --help

OR you can insert the following lines into your template:

.. code-block:: yaml

   template:
     id: pawsey
     container_dir: <value>



Fields
-------

**Required**

=============  =============  ==================================================
ID             Type           Documentation
=============  =============  ==================================================
container_dir  <class 'str'>  Location where to save and execute containers from
=============  =============  ==================================================

**Optional**

==============================  ===================================  ==========================================  ==========================================================================================
ID                              Type                                 Default                                     Documentation
==============================  ===================================  ==========================================  ==========================================================================================
intermediate_execution_dir      <class 'str'>
queues                          typing.Union[str, typing.List[str]]  workq                                       A single or list of queues that woork should be submitted to
singularity_version             <class 'str'>                        3.3.0                                       Version of singularity to load
catch_slurm_errors              <class 'bool'>                       True                                        Catch Slurm errors (like OOM or walltime)
send_job_emails                 <class 'bool'>                       True                                        (requires JanisConfiguration.notifications.email to be set) Send emails for mail types END
singularity_build_instructions  <class 'str'>                        singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
max_cores                       <class 'int'>                        28                                          Maximum number of cores a task can request
max_ram                         <class 'int'>                        128                                         Maximum amount of ram (GB) that a task can request
submission_queue                <class 'str'>                        longq                                       Queue to submit the janis 'brain' to
max_workflow_time               <class 'int'>                        5700
janis_memory_mb
==============================  ===================================  ==========================================  ==========================================================================================

