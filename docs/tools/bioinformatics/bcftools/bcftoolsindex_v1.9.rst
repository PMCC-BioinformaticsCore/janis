:orphan:

BCFTools: Index
===============================

0 contributors · 1 version

:ID: ``bcftoolsIndex``
:Python: ``janis_bioinformatics.tools.bcftools.index.versions import BcfToolsIndex_1_9``
:Versions: v1.9
:Container: michaelfranklin/bcftools:1.9
:Authors: 
:Citations: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9
:DOI: http://www.ncbi.nlm.nih.gov/pubmed/19505943
:Created: 2019-01-24
:Updated: None
:Required inputs:
   - ``vcf: CompressedVCF``
:Outputs: 
   - ``out: CompressedIndexedVCF``

Documentation
-------------

URL: `https://samtools.github.io/bcftools/bcftools.html#norm <https://samtools.github.io/bcftools/bcftools.html#norm>`_

Index bgzip compressed VCF/BCF files for random access.

------

Additional configuration (inputs)
---------------------------------

========  =================  ============================================================
name      type               documentation
========  =================  ============================================================
vcf       CompressedVCF
csi       Optional<Boolean>  (-c) generate CSI-format index for VCF/BCF files [default]
force     Optional<Boolean>  (-f) overwrite index if it already exists
minShift  Optional<Integer>  (-m) set minimal interval size for CSI indices to 2^INT [14]
tbi       Optional<Boolean>  (-t) generate TBI-format index for VCF files
threads   Optional<Integer>  sets the number of threads [0]
nrecords  Optional<Boolean>  (-n) print number of records based on existing index file
stats     Optional<Boolean>  (-s) print per contig stats based on existing index file
========  =================  ============================================================

