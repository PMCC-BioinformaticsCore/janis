:orphan:

UncompressArchive
=================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-unix>`_


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.uncompressarchive import UncompressArchive

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "uncompressarchive_step",
           UncompressArchive(
               file=None,
           )
       )
       wf.output("out", source=uncompressarchive_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for UncompressArchive:

.. code-block:: bash

   # user inputs
   janis inputs UncompressArchive > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       file: file




5. Run UncompressArchive with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       UncompressArchive





Information
------------


:ID: ``UncompressArchive``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu:latest
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>
======  ============  ===============



Additional configuration (inputs)
---------------------------------

==========  =================  ===========  ==========  =======================================================
name        type               prefix         position  documentation
==========  =================  ===========  ==========  =======================================================
file        File                                     1
stdout      Optional<Boolean>  -c                       write on standard output, keep original files unchanged
decompress  Optional<Boolean>  -d                       decompress
force       Optional<Boolean>  -f                       force overwrite of output file and compress links
keep        Optional<Boolean>  -k                       keep (don't delete) input files
list        Optional<Boolean>  -l                       list compressed file contents
noName      Optional<Boolean>  -n                       do not save or restore the original name and time stamp
name        Optional<Boolean>  -N                       save or restore the original name and time stamp
quiet       Optional<Boolean>  -q                       suppress all warnings
recursive   Optional<Boolean>  -r                       operate recursively on directories
suffix      Optional<String>   -s                       use suffix SUF on compressed files
test        Optional<Boolean>  -t                       test compressed file integrity
fast        Optional<Boolean>  -1                       compress faster
best        Optional<Boolean>  -9                       compress better
rsyncable   Optional<Boolean>  --rsyncable              Make rsync-friendly archive
==========  =================  ===========  ==========  =======================================================
