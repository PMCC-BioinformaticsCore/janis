Singularity
===========

Template ID: ``singularity``



Quickstart
-----------

Take note below how to configure the template. This quickstart only includes the fields you absolutely require. This will write a new configuration to ``~/.janis.conf``. See `configuring janis <https://janis.readthedocs.io/en/latest/references/configuration.html>`__ for more information.

.. code-block:: bash

   janis init singularity \
       --container_dir <value>
   
   # or to find out more information
   janis init singularity --help

OR you can insert the following lines into your template:

.. code-block:: yaml

   template:
     id: singularity
     container_dir: <value>



Fields
-------

**Required**

=============  ==============  ===============
ID             Type            Documentation
=============  ==============  ===============
container_dir  <class 'type'>
=============  ==============  ===============

**Optional**

=============================  =============  ==========================================  ==================================================================================
ID                             Type           Default                                     Documentation
=============================  =============  ==========================================  ==================================================================================
singularity_load_instructions                                                             Ensure singularity with this command executed in shell
container_build_instructions   <class 'str'>  singularity pull $image docker://${docker}  Instructions for building singularity, it's recommended to not touch this setting.
mail_program                   <class 'str'>                                              Mail program to pipe email to, eg: 'sendmail -t'
=============================  =============  ==========================================  ==================================================================================

