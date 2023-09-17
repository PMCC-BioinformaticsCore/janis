

Quickstart Guide
================

Summary
-------

| To translate a tool / workflow,  we use the ``translate`` functionality of `janis`.
| If we have a local CWL tool, ``tools/samtools_flagstat.cwl``, we can translate it to Nextflow with the following command:

..  code-block:: bash

   janis translate --from cwl --to nextflow tools/samtools_flagstat.cwl

| Translating a workflow has the same syntax:

..  code-block:: bash

   janis translate --from cwl --to nextflow workflows/align_sort_markdup.cwl

| Note that we only need to point janis to the main workflow file we wish to translate. All supplementary files are detected and translated automatically.

Sample Translations
-------------------

| Try out the translation functionality by following the instructions below.
| (Go to the bottom to see the full usage of the ``translate`` command)

Install Janis

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate
   pip install janis-pipelines

Obtain Sample Files

.. code-block:: bash

   git clone https://github.com/GraceAHall/sample_translation_files

Sample Tool Translations

.. code-block:: bash

   # From CWL
   janis translate --from cwl --to nextflow sample_translation_files/cwl_tool/samtools_flagstat.cwl
   janis translate --from cwl --to wdl sample_translation_files/cwl_tool/samtools_flagstat.cwl

   # From Galaxy
   janis translate --from galaxy --to cwl sample_translation_files/galaxy_tool/samtools_flagstat.xml
   janis translate --from galaxy --to nextflow sample_translation_files/galaxy_tool/samtools_flagstat.xml
   janis translate --from galaxy --to wdl sample_translation_files/galaxy_tool/samtools_flagstat.xml
   
   # From WDL
   janis translate --from wdl --to cwl sample_translation_files/wdl_tool/cutadapt.wdl
   janis translate --from wdl --to nextflow sample_translation_files/wdl_tool/cutadapt.wdl


Sample Workflow Translations

.. code-block:: bash

   # From CWL
   janis translate --from cwl --to nextflow sample_translation_files/cwl_workflow/align_sort_markdup.cwl
   janis translate --from cwl --to wdl sample_translation_files/cwl_workflow/align_sort_markdup.cwl

   # From Galaxy
   janis translate --from galaxy --to cwl sample_translation_files/galaxy_workflow/rna-seq-reads-to-counts.ga
   janis translate --from galaxy --to nextflow sample_translation_files/galaxy_workflow/rna-seq-reads-to-counts.ga
   janis translate --from galaxy --to wdl sample_translation_files/galaxy_workflow/rna-seq-reads-to-counts.ga
   
   # From WDL
   janis translate --from wdl --to cwl sample_translation_files/wdl_workflow/linear.wdl
   janis translate --from wdl --to nextflow sample_translation_files/wdl_workflow/linear.wdl

Full Usage
----------

..  raw:: html

    <div style="font-size: 13px; padding: 15px; border: 1px solid #e1e4e5ff; color: #24292e; background-color: #f8f8f8ff; font-family: SFMono-Regular,Consolas,Liberation Mono,Menlo,monospace;"><pre>
    <b>NAME</b>
        janis translate - generate translation from a source tool / workflow 

    <b>SYNOPSIS</b>
        <b>janis translate --from</b> src <b>--to</b> dest [OPTION] infile

    <b>DESCRIPTION</b>
        Ingests a source tool / workflow to generate a translation in the destination language.

        <b>infile</b> 
            Path to the source tool / main workflow file.
            For workflow translation, all supplementary files are detected and translated. 

        <b>--from</b> <u>STRING</u>
            Source language of infile.  Options: cwl | wdl | galaxy
        
        <b>--to</b> <u>STRING</u>
            Destination language to translate to.  Options: cwl | wdl | nextflow

        <b>--mode</b> <u>STRING</u>
            Translation flavour.  Options: skeleton | regular | extended.  Default: extended
                extended: full translation as close to original as possible.
                regular:  only "required" tool inputs / outputs appear during workflow translation. 
                          inferred based on workflow connections & optionality. 
                skeleton: as regular mode, except command blocks are not templated. 
                          for CWL, InputBindings are not templated for tool inputs.
        
        <b>--output-dir</b> <u>STRING</u>
            Path to the output directory for the translated files.  Default: "./translated"

    </pre></div>



|

