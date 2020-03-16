:orphan:

Parse FastQC Adaptors
===========================================

*0 contributors · 1 version*

:ID: ``ParseFastqcAdaptors``
:Python: ``janis_bioinformatics.tools.pmac.parsefastqc.v0_1_0 import ParseFastqcAdaptors``
:Versions: v0.1.0
:Container: python:3.8.1
:Authors: 
:Citations: None
:Created: 2020-01-07
:Updated: None
:Required inputs:
   - ``fastqc_datafiles: Array<File>``
:Outputs: 
   - ``adaptor_sequences: Array<String>``

Documentation
-------------

URL: *No URL to the documentation was provided*

Parse overrepresented region and lookup in Cutadapt table

Additional configuration (inputs)
---------------------------------

========================  ==============  ==========================================================================================
name                      type            prefix
========================  ==============  ==========================================================================================
fastqc_datafiles          Array<File>
cutadapt_adaptors_lookup  Optional<File>  Specifies a file which contains the list of adapter sequences which will
                                          be explicity searched against the library. The file must contain sets of named adapters in
                                          the form name[tab]sequence. Lines prefixed with a hash will be ignored.
========================  ==============  ==========================================================================================

