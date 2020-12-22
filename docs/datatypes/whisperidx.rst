
WhisperIdx
==========

A local file

Secondary Files
---------------

- ``.whisper_idx.lut_long_dir``
- ``.whisper_idx.lut_long_rc``
- ``.whisper_idx.lut_short_dir``
- ``.whisper_idx.lut_short_rc``
- ``.whisper_idx.ref_seq_desc``
- ``.whisper_idx.ref_seq_dir_pck``
- ``.whisper_idx.ref_seq_rc_pck``
- ``.whisper_idx.sa_dir``
- ``.whisper_idx.sa_rc``

.. note:: 

   For more information, visit `Secondary / Accessory Files <https://janis.readthedocs.io/en/latest/references/secondaryfiles.html>`__


Quickstart
-----------

.. code-block:: python

   from janis_bioinformatics.tools.whisper.data_types import WhisperIdx

   w = WorkflowBuilder("my_workflow")

   w.input("input_whisperidx", WhisperIdx(optional=False))
   
   # ...other workflow steps

*This page was automatically generated on 2020-12-22*.
