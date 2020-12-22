:orphan:

BCFTools: Index
===============================

``bcftoolsIndex`` · *1 contributor · 1 version*

Index bgzip compressed VCF/BCF files for random access.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bcftools.index.versions import BcfToolsIndex_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcftoolsindex_step",
           BcfToolsIndex_1_9(
               vcf=None,
           )
       )
       wf.output("out", source=bcftoolsindex_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bcftoolsIndex:

.. code-block:: bash

   # user inputs
   janis inputs bcftoolsIndex > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run bcftoolsIndex with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bcftoolsIndex





Information
------------

:ID: ``bcftoolsIndex``
:URL: `https://samtools.github.io/bcftools/bcftools.html#norm <https://samtools.github.io/bcftools/bcftools.html#norm>`_
:Versions: v1.9
:Container: biocontainers/bcftools:v1.9-1-deb_cv1
:Authors: Michael Franklin
:Citations: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9
:DOI: http://www.ncbi.nlm.nih.gov/pubmed/19505943
:Created: 2019-01-24
:Updated: 2019-01-24


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     Gzipped<VCF>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

========  =================  ===========  ==========  ============================================================
name      type               prefix         position  documentation
========  =================  ===========  ==========  ============================================================
vcf       Gzipped<VCF>                             1
csi       Optional<Boolean>  --csi                    (-c) generate CSI-format index for VCF/BCF files [default]
force     Optional<Boolean>  --force                  (-f) overwrite index if it already exists
minShift  Optional<Integer>  --min-shift              (-m) set minimal interval size for CSI indices to 2^INT [14]
tbi       Optional<Boolean>  --tbi                    (-t) generate TBI-format index for VCF files
threads   Optional<Integer>  --threads                sets the number of threads [0]
nrecords  Optional<Boolean>  --nrecords               (-n) print number of records based on existing index file
stats     Optional<Boolean>  --stats                  (-s) print per contig stats based on existing index file
========  =================  ===========  ==========  ============================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bcftoolsIndex {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       Boolean? csi
       Boolean? force
       Int? minShift
       Boolean? tbi
       Int? threads
       Boolean? nrecords
       Boolean? stats
     }
     command <<<
       set -e
       cp -f '~{vcf}' '.'
       bcftools index \
         ~{if (defined(csi) && select_first([csi])) then "--csi" else ""} \
         ~{if (defined(force) && select_first([force])) then "--force" else ""} \
         ~{if defined(minShift) then ("--min-shift " + minShift) else ''} \
         ~{if select_first([tbi, true]) then "--tbi" else ""} \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("--threads " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if (defined(nrecords) && select_first([nrecords])) then "--nrecords" else ""} \
         ~{if (defined(stats) && select_first([stats])) then "--stats" else ""} \
         '~{basename(vcf)}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = basename(vcf)
       File out_tbi = basename(vcf) + ".tbi"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'BCFTools: Index'
   doc: Index bgzip compressed VCF/BCF files for random access.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.vcf)
   - class: DockerRequirement
     dockerPull: biocontainers/bcftools:v1.9-1-deb_cv1

   inputs:
   - id: vcf
     label: vcf
     type: File
     inputBinding:
       position: 1
   - id: csi
     label: csi
     doc: (-c) generate CSI-format index for VCF/BCF files [default]
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --csi
   - id: force
     label: force
     doc: (-f) overwrite index if it already exists
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force
   - id: minShift
     label: minShift
     doc: (-m) set minimal interval size for CSI indices to 2^INT [14]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-shift
   - id: tbi
     label: tbi
     doc: (-t) generate TBI-format index for VCF files
     type: boolean
     default: true
     inputBinding:
       prefix: --tbi
   - id: threads
     label: threads
     doc: sets the number of threads [0]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads
       valueFrom: $([inputs.runtime_cpu, 1].filter(function (inner) { return inner !=
         null })[0])
   - id: nrecords
     label: nrecords
     doc: (-n) print number of records based on existing index file
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --nrecords
   - id: stats
     label: stats
     doc: (-s) print per contig stats based on existing index file
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --stats

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $(inputs.vcf.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bcftools
   - index
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bcftoolsIndex


