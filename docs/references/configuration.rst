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

   In this guide we might use dots (`.`) to refer to a nested key. For example, `notifications.email` refers to the following structure:

   .. code-block:: yaml

      notifications:
        email: <value>


Location of config
====================

By default,

- Your configuration directory is placed at `~/.janis/`.
- Your configuration path is a file called `janis.conf` in this directory, eg: `~/.janis/janis.conf`


Both `janis run / translate` allow you to provide a location to a config using the `-c / --config` parameter.

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

.. autoclass:: janis_assistant.management.envvariables.EnvVariables
   :members:


Configuration keys
===================

We use a class definition to automatically document which keys you can provide to build a Janis configuration. These keys exactly match the keys you should provide in your YAML dictionary. The type is either a Python literal, or a dictionary.

.. autoclass:: janis_assistant.management.configuration.JanisConfiguration(type)
   :members: __init__

Template
+++++++++

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationTemplate(type)
   :members: __init__


Cromwell
+++++++++

Sometimes Cromwell can be hard to configure from a simple template, so we've exposed some extra common options. Feel free to [raise an issue](https://github.com/PMCC-BioinformaticsCore/janis-assistant/issues/new) if you have questions or ideas.

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationCromwell(type)
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

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationRecipes(type)
   :members: __init__

Notifications
++++++++++++++

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationNotifications(type)
   :members: __init__

Environment
++++++++++++++

.. autoclass:: janis_assistant.management.configuration.JanisConfigurationEnvironment(type)
   :members: __init__

