:orphan:

BEDTools: intersectBed
=============================================

``bedtoolsintersectBed`` · *1 contributor · 1 version*

By far, the most common question asked of two sets of genomic features is whether or not any of the features in the two sets “overlap” with one another. This is known as feature intersection. bedtools intersect allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have fine control as to how the intersections are reported. bedtools intersect works with both BED/GFF/VCF and BAM files as input.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bedtools.intersectbed.versions import BedToolsIntersectBed_2_29_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bedtoolsintersectbed_step",
           BedToolsIntersectBed_2_29_2(
               inputABam=None,
               inputBBed=None,
           )
       )
       wf.output("out", source=bedtoolsintersectbed_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bedtoolsintersectBed:

.. code-block:: bash

   # user inputs
   janis inputs bedtoolsintersectBed > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputABam: inputABam.bam
       inputBBed:
       - inputBBed_0.bed
       - inputBBed_1.bed




5. Run bedtoolsintersectBed with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bedtoolsintersectBed





Information
------------

:ID: ``bedtoolsintersectBed``
:URL: `https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html <https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html>`_
:Versions: v2.29.2
:Container: quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-02-20
:Updated: 2020-02-26


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<BAM>
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

=====================  =================  ============  ==========  ===========================================================================================================================================================================================================================================================================================================================================
name                   type               prefix        position    documentation
=====================  =================  ============  ==========  ===========================================================================================================================================================================================================================================================================================================================================
inputABam              BAM                -a                        input file a: only bam is supported at the moment
inputBBed              Array<bed>         -b                        input file b: only bed is supported at the moment. May be followed with multiple databases and/or  wildcard (*) character(s).
writeOriginalA         Optional<Boolean>  -wa                       Write the original entry in A for each overlap.
writeOriginalB         Optional<Boolean>  -wb                       Write the original entry in B for each overlap. - Useful for knowing _what_ A overlaps. Restricted by -f  and -r.
leftOuterJoin          Optional<Boolean>  -loj                      Perform a 'left outer join'. That is, for each feature in A report each overlap with B.  If no overlaps are found, report a NULL feature for B.
writeOriginalAB        Optional<Boolean>  -wo                       Write the original A and B entries plus the number of base pairs of overlap between the two features. - Overlaps restricted by -f and -r. Only A features with overlap are reported.
writeABBase            Optional<Boolean>  -wao                      Write the original A and B entries plus the number of base pairs of overlap between the two features. - Overlapping features restricted by -f and -r. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
modeu                  Optional<Boolean>  -u                        Write the original A entry _once_ if _any_ overlaps found in B. - In other words, just report the fact >=1 hit was found. - Overlaps restricted by -f and -r.
modec                  Optional<Boolean>  -c                        For each entry in A, report the number of overlaps with B. - Reports 0 for A entries that have no overlap with B. - Overlaps restricted by -f, -F, -r, and -s.
modeC                  Optional<Boolean>  -C                        -C	For each entry in A, separately report the number of - overlaps with each B file on a distinct line. - Reports 0 for A entries that have no overlap with B. - Overlaps restricted by -f, -F, -r, and -s.
modev                  Optional<Boolean>  -v                        Only report those entries in A that have _no overlaps_ with B. - Similar to 'grep -v' (an homage).
strandedness           Optional<Boolean>  -s                        Require same strandedness.  That is, only report hits in B that overlap A on the _same_ strand. - By default, overlaps are reported without respect to strand.
differentStrandedness  Optional<Boolean>  -S                        Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand. - By default, overlaps are reported without respect to strand.
fractionA              Optional<Float>    -f                        Minimum overlap required as a fraction of A. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
fractionB              Optional<Float>    -F                        Minimum overlap required as a fraction of B. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
reciprocalFraction     Optional<Boolean>  -r                        Require that the fraction overlap be reciprocal for A AND B. - In other words, if -f is 0.90 and -r is used, this requires that B overlap 90% of A and A _also_ overlaps 90% of B.
minFraction            Optional<Boolean>  -r                        Require that the minimum fraction be satisfied for A OR B. - In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of  B is covered. Without -e, both fractions would have to be satisfied.
split                  Optional<Boolean>  -split                    Treat 'split' BAM or BED12 entries as distinct BED intervals.
genome                 Optional<File>     -g                        Provide a genome file to enforce consistent chromosome sort order across input files. Only applies when used with -sorted option.
noNameCheck            Optional<Boolean>  -nonamecheck              For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. 'chr1' vs 'chr01'.
sorted                 Optional<Boolean>  -sorted                   Use the 'chromsweep' algorithm for sorted (-k1,1 -k2,2n) input.
sortOut                Optional<Boolean>  -sortout                  When using multiple databases, sort the output DB hits for each record.
header                 Optional<Boolean>  -header                   Print the header from the A file prior to results.
noBuf                  Optional<Boolean>  -nobuf                    Disable buffered output. Using this option will cause each line of output to be printed as it is generated, rather than saved in a buffer. This will make printing large output files noticeably slower, but can be useful in conjunction with other software tools and scripts that need to process one line of bedtools output at a time.
bufMem                 Optional<Integer>  -iobuf                    Specify amount of memory to use for input buffer. Takes an integer argument. Optional suffixes K/M/G supported. Note: currently has no effect with compressed files.
=====================  =================  ============  ==========  ===========================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bedtoolsintersectBed {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Boolean? writeOriginalA
       Boolean? writeOriginalB
       Boolean? leftOuterJoin
       Boolean? writeOriginalAB
       Boolean? writeABBase
       Boolean? modeu
       Boolean? modec
       Boolean? modeC
       Boolean? modev
       Boolean? strandedness
       Boolean? differentStrandedness
       Float? fractionA
       Float? fractionB
       Boolean? reciprocalFraction
       Boolean? minFraction
       Boolean? split
       File? genome
       Boolean? noNameCheck
       Boolean? sorted
       Boolean? sortOut
       Boolean? header
       Boolean? noBuf
       Int? bufMem
       File inputABam
       Array[File] inputBBed
     }
     command <<<
       set -e
       intersectBed \
         ~{if (defined(writeOriginalA) && select_first([writeOriginalA])) then "-wa" else ""} \
         ~{if (defined(writeOriginalB) && select_first([writeOriginalB])) then "-wb" else ""} \
         ~{if (defined(leftOuterJoin) && select_first([leftOuterJoin])) then "-loj" else ""} \
         ~{if (defined(writeOriginalAB) && select_first([writeOriginalAB])) then "-wo" else ""} \
         ~{if (defined(writeABBase) && select_first([writeABBase])) then "-wao" else ""} \
         ~{if (defined(modeu) && select_first([modeu])) then "-u" else ""} \
         ~{if (defined(modec) && select_first([modec])) then "-c" else ""} \
         ~{if (defined(modeC) && select_first([modeC])) then "-C" else ""} \
         ~{if (defined(modev) && select_first([modev])) then "-v" else ""} \
         ~{if (defined(strandedness) && select_first([strandedness])) then "-s" else ""} \
         ~{if (defined(differentStrandedness) && select_first([differentStrandedness])) then "-S" else ""} \
         ~{if defined(fractionA) then ("-f " + fractionA) else ''} \
         ~{if defined(fractionB) then ("-F " + fractionB) else ''} \
         ~{if (defined(reciprocalFraction) && select_first([reciprocalFraction])) then "-r" else ""} \
         ~{if (defined(minFraction) && select_first([minFraction])) then "-r" else ""} \
         ~{if (defined(split) && select_first([split])) then "-split" else ""} \
         ~{if defined(genome) then ("-g '" + genome + "'") else ""} \
         ~{if (defined(noNameCheck) && select_first([noNameCheck])) then "-nonamecheck" else ""} \
         ~{if (defined(sorted) && select_first([sorted])) then "-sorted" else ""} \
         ~{if (defined(sortOut) && select_first([sortOut])) then "-sortout" else ""} \
         ~{if (defined(header) && select_first([header])) then "-header" else ""} \
         ~{if (defined(noBuf) && select_first([noBuf])) then "-nobuf" else ""} \
         ~{if defined(bufMem) then ("-iobuf " + bufMem) else ''} \
         -a '~{inputABam}' \
         ~{if length(inputBBed) > 0 then "-b '" + sep("' '", inputBBed) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'BEDTools: intersectBed'
   doc: |-
     By far, the most common question asked of two sets of genomic features is whether or not any of the features in the two sets “overlap” with one another. This is known as feature intersection. bedtools intersect allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have fine control as to how the intersections are reported. bedtools intersect works with both BED/GFF/VCF and BAM files as input.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0

   inputs:
   - id: writeOriginalA
     label: writeOriginalA
     doc: Write the original entry in A for each overlap.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -wa
   - id: writeOriginalB
     label: writeOriginalB
     doc: |-
       Write the original entry in B for each overlap. - Useful for knowing _what_ A overlaps. Restricted by -f  and -r.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -wb
   - id: leftOuterJoin
     label: leftOuterJoin
     doc: |-
       Perform a 'left outer join'. That is, for each feature in A report each overlap with B.  If no overlaps are found, report a NULL feature for B.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -loj
   - id: writeOriginalAB
     label: writeOriginalAB
     doc: |-
       Write the original A and B entries plus the number of base pairs of overlap between the two features. - Overlaps restricted by -f and -r. Only A features with overlap are reported.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -wo
   - id: writeABBase
     label: writeABBase
     doc: |-
       Write the original A and B entries plus the number of base pairs of overlap between the two features. - Overlapping features restricted by -f and -r. However, A features w/o overlap are also reported with a NULL B feature and overlap = 0.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -wao
   - id: modeu
     label: modeu
     doc: |-
       Write the original A entry _once_ if _any_ overlaps found in B. - In other words, just report the fact >=1 hit was found. - Overlaps restricted by -f and -r.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -u
   - id: modec
     label: modec
     doc: |-
       For each entry in A, report the number of overlaps with B. - Reports 0 for A entries that have no overlap with B. - Overlaps restricted by -f, -F, -r, and -s.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -c
   - id: modeC
     label: modeC
     doc: |-
       -C	For each entry in A, separately report the number of - overlaps with each B file on a distinct line. - Reports 0 for A entries that have no overlap with B. - Overlaps restricted by -f, -F, -r, and -s.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -C
   - id: modev
     label: modev
     doc: |-
       Only report those entries in A that have _no overlaps_ with B. - Similar to 'grep -v' (an homage).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -v
   - id: strandedness
     label: strandedness
     doc: |-
       Require same strandedness.  That is, only report hits in B that overlap A on the _same_ strand. - By default, overlaps are reported without respect to strand.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -s
   - id: differentStrandedness
     label: differentStrandedness
     doc: |-
       Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand. - By default, overlaps are reported without respect to strand.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -S
   - id: fractionA
     label: fractionA
     doc: |-
       Minimum overlap required as a fraction of A. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -f
   - id: fractionB
     label: fractionB
     doc: |-
       Minimum overlap required as a fraction of B. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -F
   - id: reciprocalFraction
     label: reciprocalFraction
     doc: |-
       Require that the fraction overlap be reciprocal for A AND B. - In other words, if -f is 0.90 and -r is used, this requires that B overlap 90% of A and A _also_ overlaps 90% of B.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -r
   - id: minFraction
     label: minFraction
     doc: |-
       Require that the minimum fraction be satisfied for A OR B. - In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of  B is covered. Without -e, both fractions would have to be satisfied.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -r
   - id: split
     label: split
     doc: Treat 'split' BAM or BED12 entries as distinct BED intervals.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -split
   - id: genome
     label: genome
     doc: |-
       Provide a genome file to enforce consistent chromosome sort order across input files. Only applies when used with -sorted option.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -g
   - id: noNameCheck
     label: noNameCheck
     doc: |-
       For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. 'chr1' vs 'chr01'.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -nonamecheck
   - id: sorted
     label: sorted
     doc: Use the 'chromsweep' algorithm for sorted (-k1,1 -k2,2n) input.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -sorted
   - id: sortOut
     label: sortOut
     doc: When using multiple databases, sort the output DB hits for each record.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -sortout
   - id: header
     label: header
     doc: Print the header from the A file prior to results.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -header
   - id: noBuf
     label: noBuf
     doc: |-
       Disable buffered output. Using this option will cause each line of output to be printed as it is generated, rather than saved in a buffer. This will make printing large output files noticeably slower, but can be useful in conjunction with other software tools and scripts that need to process one line of bedtools output at a time.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -nobuf
   - id: bufMem
     label: bufMem
     doc: |-
       Specify amount of memory to use for input buffer. Takes an integer argument. Optional suffixes K/M/G supported. Note: currently has no effect with compressed files.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -iobuf
   - id: inputABam
     label: inputABam
     doc: 'input file a: only bam is supported at the moment'
     type: File
     inputBinding:
       prefix: -a
   - id: inputBBed
     label: inputBBed
     doc: |-
       input file b: only bed is supported at the moment. May be followed with multiple databases and/or  wildcard (*) character(s). 
     type:
       type: array
       items: File
     inputBinding:
       prefix: -b

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - intersectBed
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bedtoolsintersectBed


