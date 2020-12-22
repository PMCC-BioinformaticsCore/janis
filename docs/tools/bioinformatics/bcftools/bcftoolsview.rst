:orphan:

BCFTools: View
=============================

``bcftoolsview`` · *1 contributor · 2 versions*

________________________________
 
        View, subset and filter VCF or BCF files by position and filtering expression
        Convert between VCF and BCF. Former bcftools subset.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bcftools.view.versions import BcfToolsView_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcftoolsview_step",
           BcfToolsView_1_9(
               file=None,
           )
       )
       wf.output("out", source=bcftoolsview_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bcftoolsview:

.. code-block:: bash

   # user inputs
   janis inputs bcftoolsview > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       file: file.vcf.gz




5. Run bcftoolsview with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bcftoolsview





Information
------------

:ID: ``bcftoolsview``
:URL: `https://samtools.github.io/bcftools/bcftools.html#view <https://samtools.github.io/bcftools/bcftools.html#view>`_
:Versions: v1.9, v1.5
:Container: biocontainers/bcftools:v1.9-1-deb_cv1
:Authors: Michael Franklin
:Citations: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9
:DOI: http://www.ncbi.nlm.nih.gov/pubmed/19505943
:Created: 2019-01-24
:Updated: 2019-01-24


Outputs
-----------

======  ====================  ===============
name    type                  documentation
======  ====================  ===============
out     stdout<Gzipped<VCF>>
======  ====================  ===============


Additional configuration (inputs)
---------------------------------

================  =======================  ===================  ==========  ==============================================================================================================================================================================
name              type                     prefix                 position  documentation
================  =======================  ===================  ==========  ==============================================================================================================================================================================
file              Gzipped<VCF>                                           2
dropGenotypes     Optional<Boolean>        --drop-genotypes              1  (-G) drop individual genotype information (after subsetting if -s option set)
headerOnly        Optional<Boolean>        --header-only                 1  (-h) print the header only
noHeader          Optional<Boolean>        --no-header                   1  (-H) suppress the header in VCF output
compressionLevel  Optional<Integer>        --compression-level           1  (-l) compression level: 0 uncompressed, 1 best speed, 9 best compression [-1]
noVersion         Optional<Boolean>        --no-version                  1  do not append version and command line to the header
regions           Optional<String>         --regions                     1  (-r) restrict to comma-separated list of regions
regionsFile       Optional<File>           --regions-file                1  (-R) restrict to regions listed in a file
targets           Optional<String>         --targets                     1  (-t) similar to -r but streams rather than index-jumps. Exclude regions with '^' prefix
targetsFile       Optional<File>           --targets-file                1  (-T) similar to -R but streams rather than index-jumps. Exclude regions with '^' prefix
threads           Optional<Integer>        --threads                     1  number of extra output compression threads [0]
trimAltAlleles    Optional<Boolean>        --trim-alt-alleles            1  (-a) trim alternate alleles not seen in the subset
noUpdate          Optional<Boolean>        --no-update                   1  (-I) do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)
samples           Optional<Array<String>>  --samples                     1  (-s) comma separated list of samples to include (or exclude with '^' prefix)
samplesFile       Optional<File>           --samples-file                1  (-S) file of samples to include (or exclude with '^' prefix)
forceSamples      Optional<Boolean>        --force-samples               1  only warn about unknown subset samples
minAc             Optional<Integer>        --min-ac                      1  (-c) minimum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
maxAc             Optional<Integer>        --max-ac                      1  (-C) maximum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
applyFilters      Optional<Array<String>>  --apply-filters               1  (-f) require at least one of the listed FILTER strings (e.g. 'PASS,.'')
genotype          Optional<String>         --genotype                    1  (-g) [<hom|het|miss>] require one or more hom/het/missing genotype or, if prefixed with '^', exclude sites with hom/het/missing genotypes
include           Optional<String>         --include                     1  (-i) select sites for which the expression is true (see man page for details)
exclude           Optional<String>         --exclude                     1  (-e) exclude sites for which the expression is true (see man page for details)
known             Optional<Boolean>        --known                       1  (-k) select known sites only (ID is not/is '.')
novel             Optional<Boolean>        --novel                       1  (-n) select novel sites only (ID is not/is '.')
minAlleles        Optional<Integer>        --min-alleles                 1  (-m) minimum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
maxAlleles        Optional<Integer>        --max-alleles                 1  (-M) maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
phased            Optional<Boolean>        --phased                      1  (-p) select sites where all samples are phased
excludePhased     Optional<Boolean>        --exclude-phased              1  (-P) exclude sites where all samples are phased
minAf             Optional<Float>          --min-af                      1  (-q) minimum frequency for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
maxAf             Optional<Float>          --max-af                      1  (-Q) maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
uncalled          Optional<Boolean>        --uncalled                    1  (-u) select sites without a called genotype
excludeUncalled   Optional<Boolean>        --exclude-uncalled            1  (-U) exclude sites without a called genotype
types             Optional<Array<String>>  --types                       1  (-v) select comma-separated list of variant types: snps,indels,mnps,other [null]
excludeTypes      Optional<Array<String>>  --exclude-types               1  (-V) exclude comma-separated list of variant types: snps,indels,mnps,other [null]
private           Optional<Boolean>        --private                     1  (-x) select sites where the non-reference alleles are exclusive (private) to the subset samples
excludePrivate    Optional<Boolean>        --exclude-private             1  (-X) exclude sites where the non-reference alleles are exclusive (private) to the subset samples
================  =======================  ===================  ==========  ==============================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bcftoolsview {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File file
       Boolean? dropGenotypes
       Boolean? headerOnly
       Boolean? noHeader
       Int? compressionLevel
       Boolean? noVersion
       String? regions
       File? regionsFile
       String? targets
       File? targetsFile
       Int? threads
       Boolean? trimAltAlleles
       Boolean? noUpdate
       Array[String]? samples
       File? samplesFile
       Boolean? forceSamples
       Int? minAc
       Int? maxAc
       Array[String]? applyFilters
       String? genotype
       String? include
       String? exclude
       Boolean? known
       Boolean? novel
       Int? minAlleles
       Int? maxAlleles
       Boolean? phased
       Boolean? excludePhased
       Float? minAf
       Float? maxAf
       Boolean? uncalled
       Boolean? excludeUncalled
       Array[String]? types
       Array[String]? excludeTypes
       Boolean? private
       Boolean? excludePrivate
     }
     command <<<
       set -e
       bcftools view \
         ~{if (defined(dropGenotypes) && select_first([dropGenotypes])) then "--drop-genotypes" else ""} \
         ~{if (defined(headerOnly) && select_first([headerOnly])) then "--header-only" else ""} \
         ~{if (defined(noHeader) && select_first([noHeader])) then "--no-header" else ""} \
         ~{if defined(compressionLevel) then ("--compression-level " + compressionLevel) else ''} \
         ~{if (defined(noVersion) && select_first([noVersion])) then "--no-version" else ""} \
         ~{if defined(regions) then ("--regions '" + regions + "'") else ""} \
         ~{if defined(regionsFile) then ("--regions-file '" + regionsFile + "'") else ""} \
         ~{if defined(targets) then ("--targets '" + targets + "'") else ""} \
         ~{if defined(targetsFile) then ("--targets-file '" + targetsFile + "'") else ""} \
         ~{if defined(threads) then ("--threads " + threads) else ''} \
         ~{if (defined(trimAltAlleles) && select_first([trimAltAlleles])) then "--trim-alt-alleles" else ""} \
         ~{if (defined(noUpdate) && select_first([noUpdate])) then "--no-update" else ""} \
         ~{if (defined(samples) && length(select_first([samples])) > 0) then "--samples '" + sep("' '", select_first([samples])) + "'" else ""} \
         ~{if defined(samplesFile) then ("--samples-file '" + samplesFile + "'") else ""} \
         ~{if (defined(forceSamples) && select_first([forceSamples])) then "--force-samples" else ""} \
         ~{if defined(minAc) then ("--min-ac " + minAc) else ''} \
         ~{if defined(maxAc) then ("--max-ac " + maxAc) else ''} \
         ~{if (defined(applyFilters) && length(select_first([applyFilters])) > 0) then "--apply-filters '" + sep("' '", select_first([applyFilters])) + "'" else ""} \
         ~{if defined(genotype) then ("--genotype '" + genotype + "'") else ""} \
         ~{if defined(include) then ("--include '" + include + "'") else ""} \
         ~{if defined(exclude) then ("--exclude '" + exclude + "'") else ""} \
         ~{if (defined(known) && select_first([known])) then "--known" else ""} \
         ~{if (defined(novel) && select_first([novel])) then "--novel" else ""} \
         ~{if defined(minAlleles) then ("--min-alleles " + minAlleles) else ''} \
         ~{if defined(maxAlleles) then ("--max-alleles " + maxAlleles) else ''} \
         ~{if (defined(phased) && select_first([phased])) then "--phased" else ""} \
         ~{if (defined(excludePhased) && select_first([excludePhased])) then "--exclude-phased" else ""} \
         ~{if defined(minAf) then ("--min-af " + minAf) else ''} \
         ~{if defined(maxAf) then ("--max-af " + maxAf) else ''} \
         ~{if (defined(uncalled) && select_first([uncalled])) then "--uncalled" else ""} \
         ~{if (defined(excludeUncalled) && select_first([excludeUncalled])) then "--exclude-uncalled" else ""} \
         ~{if (defined(types) && length(select_first([types])) > 0) then "--types '" + sep("' '", select_first([types])) + "'" else ""} \
         ~{if (defined(excludeTypes) && length(select_first([excludeTypes])) > 0) then "--exclude-types '" + sep("' '", select_first([excludeTypes])) + "'" else ""} \
         ~{if (defined(private) && select_first([private])) then "--private" else ""} \
         ~{if (defined(excludePrivate) && select_first([excludePrivate])) then "--exclude-private" else ""} \
         --output-type 'z' \
         '~{file}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
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
   label: 'BCFTools: View'
   doc: |-
     ________________________________
   
             View, subset and filter VCF or BCF files by position and filtering expression
             Convert between VCF and BCF. Former bcftools subset.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/bcftools:v1.9-1-deb_cv1

   inputs:
   - id: file
     label: file
     type: File
     inputBinding:
       position: 2
   - id: dropGenotypes
     label: dropGenotypes
     doc: (-G) drop individual genotype information (after subsetting if -s option set)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --drop-genotypes
       position: 1
   - id: headerOnly
     label: headerOnly
     doc: (-h) print the header only
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --header-only
       position: 1
   - id: noHeader
     label: noHeader
     doc: (-H) suppress the header in VCF output
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-header
       position: 1
   - id: compressionLevel
     label: compressionLevel
     doc: '(-l) compression level: 0 uncompressed, 1 best speed, 9 best compression [-1]'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --compression-level
       position: 1
   - id: noVersion
     label: noVersion
     doc: do not append version and command line to the header
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-version
       position: 1
   - id: regions
     label: regions
     doc: (-r) restrict to comma-separated list of regions
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --regions
       position: 1
   - id: regionsFile
     label: regionsFile
     doc: (-R) restrict to regions listed in a file
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --regions-file
       position: 1
   - id: targets
     label: targets
     doc: |-
       (-t) similar to -r but streams rather than index-jumps. Exclude regions with '^' prefix
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --targets
       position: 1
   - id: targetsFile
     label: targetsFile
     doc: |-
       (-T) similar to -R but streams rather than index-jumps. Exclude regions with '^' prefix
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --targets-file
       position: 1
   - id: threads
     label: threads
     doc: number of extra output compression threads [0]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads
       position: 1
   - id: trimAltAlleles
     label: trimAltAlleles
     doc: (-a) trim alternate alleles not seen in the subset
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --trim-alt-alleles
       position: 1
   - id: noUpdate
     label: noUpdate
     doc: |-
       (-I) do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-update
       position: 1
   - id: samples
     label: samples
     doc: (-s) comma separated list of samples to include (or exclude with '^' prefix)
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --samples
       position: 1
   - id: samplesFile
     label: samplesFile
     doc: (-S) file of samples to include (or exclude with '^' prefix)
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --samples-file
       position: 1
   - id: forceSamples
     label: forceSamples
     doc: only warn about unknown subset samples
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force-samples
       position: 1
   - id: minAc
     label: minAc
     doc: |-
       (-c) minimum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-ac
       position: 1
   - id: maxAc
     label: maxAc
     doc: |-
       (-C) maximum count for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-ac
       position: 1
   - id: applyFilters
     label: applyFilters
     doc: (-f) require at least one of the listed FILTER strings (e.g. 'PASS,.'')
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --apply-filters
       position: 1
   - id: genotype
     label: genotype
     doc: |-
       (-g) [<hom|het|miss>] require one or more hom/het/missing genotype or, if prefixed with '^', exclude sites with hom/het/missing genotypes
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genotype
       position: 1
   - id: include
     label: include
     doc: (-i) select sites for which the expression is true (see man page for details)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --include
       position: 1
   - id: exclude
     label: exclude
     doc: (-e) exclude sites for which the expression is true (see man page for details)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --exclude
       position: 1
   - id: known
     label: known
     doc: (-k) select known sites only (ID is not/is '.')
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --known
       position: 1
   - id: novel
     label: novel
     doc: (-n) select novel sites only (ID is not/is '.')
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --novel
       position: 1
   - id: minAlleles
     label: minAlleles
     doc: |-
       (-m) minimum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-alleles
       position: 1
   - id: maxAlleles
     label: maxAlleles
     doc: |-
       (-M) maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-alleles
       position: 1
   - id: phased
     label: phased
     doc: (-p) select sites where all samples are phased
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --phased
       position: 1
   - id: excludePhased
     label: excludePhased
     doc: (-P) exclude sites where all samples are phased
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude-phased
       position: 1
   - id: minAf
     label: minAf
     doc: |-
       (-q) minimum frequency for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --min-af
       position: 1
   - id: maxAf
     label: maxAf
     doc: |-
       (-Q) maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --max-af
       position: 1
   - id: uncalled
     label: uncalled
     doc: (-u) select sites without a called genotype
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --uncalled
       position: 1
   - id: excludeUncalled
     label: excludeUncalled
     doc: (-U) exclude sites without a called genotype
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude-uncalled
       position: 1
   - id: types
     label: types
     doc: '(-v) select comma-separated list of variant types: snps,indels,mnps,other
       [null]'
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --types
       position: 1
   - id: excludeTypes
     label: excludeTypes
     doc: |-
       (-V) exclude comma-separated list of variant types: snps,indels,mnps,other [null]
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --exclude-types
       position: 1
   - id: private
     label: private
     doc: |-
       (-x) select sites where the non-reference alleles are exclusive (private) to the subset samples
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --private
       position: 1
   - id: excludePrivate
     label: excludePrivate
     doc: |-
       (-X) exclude sites where the non-reference alleles are exclusive (private) to the subset samples
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude-private
       position: 1

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bcftools
   - view
   arguments:
   - prefix: --output-type
     position: 1
     valueFrom: z

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bcftoolsview


