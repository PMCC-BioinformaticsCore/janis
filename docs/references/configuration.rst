Configuring Janis Assistant
##############################

.. warning:: Configuring Janis had backwards incompatible changes in v0.12.0 with regards to specific keys.

When we talk about configuring Janis, we're really talking about how to configure the `Janis assistant <https://github.com/PMCC-BioinformaticsCore/janis-assistant>`_ and hence Cromwell / CWLTool to interact with batch systems (eg: Slurm / PBS Torque) and container environments (Docker / Singularity).

Janis has built in templates for the following compute environments:

- Local (Docker or Singularity)
    - The only environment compatible with CWLTool.
- Slurm (Singularity only)
- PBS / Torque (Singularity only)

In an extension of Janis (`janis-templates <https://github.com/PMCC-BioinformaticsCore/janis-templates>`_), we've produced a number of location specific templates, usually with sensible defaults.

See this list for a full list of templates: `https://janis.readthedocs.io/en/latest/templates/index.html <https://janis.readthedocs.io/en/latest/templates/index.html>`_.

Syntax
==========

The config should be in YAML, and can contain nested dictionaries.

.. note::

   In this guide we might use dots (``.``) to refer to a nested key. For example, ``notifications.email`` refers to the following structure:

   .. code-block:: yaml

      notifications:
        email: <value>


Location of config
====================

By default,

- Your configuration directory is placed at ``~/.janis/``.
- Your configuration path is a file called ``janis.conf`` in this directory, eg: ``~/.janis/janis.conf``


Both ``janis run / translate`` allow you to provide a location to a config using the ``-c / --config`` parameter.

In addition, you can also configure the following environment variables:


- ``JANIS_CONFIGPATH`` - simply the path to your config file
- ``JANIS_CONFIGDIR`` - The configuration path is determined by `$JANIS_CONFIGDIR/janis.conf`


Initialising a template
========================

```bash
janis init <template> [...options]
```

By default this will place our config at `~/.janis/janis.conf`. If there's already a config there, it will NOT override it.

 - ``-o`` / ``--output``: output path to write to, default (`~/.janis/janis.conf`).
 - ``-f`` / ``--force``: Overwrite the config if one exists at the output path.
 - ``--stdout``: Write the output to stdout.

Environment variables
======================

The following environment variables allow you to customise the behaviour of Janis, without writing a configuration file. This is a convenient way to automatically configure Janis through an HPCs module system.

All the environment variables start with ``JANIS_``, and is the value of enums declared below (and not the key name), for example:

.. code-block:: bash

   export JANIS_CROMWELLJAR='/path/to/cromwell.jar'

.. autoclass:: janis_assistant.management.envvariables.EnvVariables
   :members:


Configuration keys
===================

We use a class definition to automatically document which keys you can provide to build a Janis configuration. These keys exactly match the keys you should provide in your YAML dictionary. The type is either a Python literal, or a dictionary.

For example, the class definition below corresponds to the following (partial) YAML configuration:

.. code-block:: yaml

   # EXAMPLE CONFIGURATION ONLY
   engine: cromwell
   run_in_background: false
   call_caching_enabled: true
   cromwell:
       jar: /Users/franklinmichael/broad/cromwell-53.1.jar
       call_caching_method: fingerprint

.. autoclass:: janis_assistant.management.configuration.JanisConfiguration(type)
   :members: __init__

Template
+++++++++

Janis templates are a convenient way to handle configuring Janis, Cromwell and CWLTool for special environments. A number of templates are prebuilt into Janis, such as ``slurm_singularity``, ``slurm_pbs``, and number of additional templates for specific HPCs (like Peter Mac, Spartan at UoM) are available, and documented in the:

- `List of templates <https://janis.readthedocs.io/en/latest/templates/index.html>`_ page.

You could use a template like the following:

.. code-block:: yaml

   # rest of janis configuration
   template:
     id: slurm_singularity
     # arguments for template 'slurm_singularity', like 'container_dir'.
     container_dir: /shared/path/to/containerdir/


.. autoclass:: janis_assistant.management.configuration.JanisConfigurationTemplate(type)
   :members: __init__


Cromwell
+++++++++

Sometimes Cromwell can be hard to configure from a simple template, so we've exposed some extra common options. Feel free to [raise an issue](https://github.com/PMCC-BioinformaticsCore/janis-assistant/issues/new) if you have questions or ideas.

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationCromwell(type)
   :members: __init__

.. autoclass:: janis_assistant.management.configuration.MySqlInstanceConfig(type)
   :members: __init__


Existing Cromwell instance
............................

In addition to the previous ``cromwell.url`` method, you can also manage this through the command line.
To configure Janis to submit to an existing Cromwell instance (eg: one already set up for the cloud), the CLI has a mechanism for setting the ``cromwell_url``:

.. code-block:: bash

   urlwithport="127.0.0.1:8000"
   janis run --engine cromwell --cromwell-url $urlwithport hello

OR

**`~/janis.conf`**

.. code-block:: yaml

   engine: cromwell
   cromwell:
     url: 127.0.0.1:8000

Overriding Cromwell JAR
.........................

In additional to the previous ``cromwell.jar``, you can set the location of the Cromwell JAR through the environment variable ``JANIS_CROMWELLJAR``:

.. code-block:: bash

   export JANIS_CROMWELLJAR=/path/to/cromwell.jar


Recipes
++++++++++

Often between a number of different workflows, you have a set of inputs that you want to apply multiple times. For that, we have `recipes`.

.. note:: Recipes only match on the input name. Your input names _MUST_ be consistent through your pipelines for this concept to be useful.

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationRecipes(type)
   :members: __init__

For example, everytime I run the `WGSGermlineGATK <https://janis.readthedocs.io/en/latest/pipelines/wgsgermlinegatk.html>`_ pipeline with ``hg38``, I know I want to provide the same reference files. There are a few ways to configure this:

- Recipes: A dictionary of input values, keyed by the recipe name.
- Paths: a list of ``*.yaml`` files, where each path contains a dictionary of input values, keyed by the recipe name, similar to the previous recipes name.
- Directories: a directory of ``*.yaml`` files, where the ``*`` is the recipe name.

The examples below, encode the following information. When we use the ``hg38`` recipe, we want to provide an input value for reference as ``/path/to/hg38/reference.fasta``, and input value for ``type`` as ``hg38``. Similar for a second recipe for ``hg19``.

.. code-block::
   hg38:
     reference: /path/to/hg38/reference.fasta
     type: hg38
   hg19:
     reference: /path/to/hg19/reference.fasta
     type: hg19


Recipes dictionary
....................

You can specify this recipe directly in your ``janis.conf``:

.. code-block:: yaml

   recipes:
     recipes:
       hg38:
         reference: /path/to/hg38/reference.fasta
         type: hg38
       hg19:
         reference: /path/to/hg19/reference.fasta
         type: hg19

Recipes Paths
..............

Or you could create a ``myrecipes.yaml`` with the contents:

.. code-block:: yaml

   hg38:
     reference: /path/to/hg38/reference.fasta
     type: hg38
   hg19:
     reference: /path/to/hg19/reference.fasta
     type: hg19

And then instruct Janis to use this file in two ways:

1. In your ``janis.conf`` with:

  .. code-block:: yaml

     recipes:
       paths:
       - /path/to/myrecipes.yaml

2. OR, you can export the comma-separated environment variable:

  .. code-block:: bash

     export JANIS_RECIPEPATHS="/path/to/myrecipes.yaml,/path/to/myrecipes2.yaml"

Recipe Directories
....................

Create two files in a directory"

1. ``hg38.yaml``:

    .. code-block:: yaml

       reference: /path/to/hg38/reference.fasta
       type: hg38

2. ``hg19.yaml``:

    .. code-block:: yaml

       reference: /path/to/hg19/reference.fasta
       type: hg19


And similar to the paths, you can specify this directory in two ways:

1. In your ``janis.conf`` with:

  .. code-block:: yaml

     recipes:
       directories:
       # /path/to/recipes has two files, hg38.yaml | hg19.yaml
       - /path/to/recipes/

2. OR, you can export the comma-separated environment variable:

  .. code-block:: bash

     export JANIS_RECIPEPATHS="/path/to/recipes/,/path/to/recipes2/"



Notifications
++++++++++++++

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationNotifications(type)
   :members: __init__

Environment
++++++++++++++

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationEnvironment(type)
   :members: __init__

