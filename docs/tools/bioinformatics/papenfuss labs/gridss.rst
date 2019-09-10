
.. include:: gridss_v2.5.1-dev
.. include:: gridss_v2.4.0
.. include:: gridss_v2.2.3

Gridss
===============

Description
-------------

Tool identifier: ``gridss``

Tool path: ``janis_bioinformatics.tools.pappenfuss.gridss.gridss import Gridss_2_5_1``

Version: v2.5.1-dev

Container: ``michaelfranklin/gridss:2.5.1-dev2``

Versions
*********

- v2.5.1-dev (current)
- `v2.4.0 <gridss_v2.4.0.html>`_
- `v2.2.3 <gridss_v2.2.3.html>`_

Documentation
-------------

URL
******
`https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation <https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation>`_

Tool documentation
******************
GRIDSS: the Genomic Rearrangement IDentification Software Suite

GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements.
GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina
sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide
break-end assembly, split read, and read pair evidence.

GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design,
any step (such as split read identification) can be replaced by another implementation that also outputs
using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant
pipeline designed for interoperability with other tools.

If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback
from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.


Outputs
-------
========  ======  ===============
name      type    documentation
========  ======  ===============
out       VCF
assembly  BAM
========  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ===========  ==========  ===============
name       type           prefix         position  documentation
=========  =============  ===========  ==========  ===============
bams       Array<BAM>                          10
reference  FastaWithDict  --reference           1
=========  =============  ===========  ==========  ===============

Optional inputs
***************

================  ==================  ===========  ==========  ===============
name              type                prefix         position  documentation
================  ==================  ===========  ==========  ===============
outputFilename    Optional<Filename>  --output              2
assemblyFilename  Optional<Filename>  --assembly            3
threads           Optional<Integer>   --threads
blacklist         Optional<bed>       --blacklist           4
================  ==================  ===========  ==========  ===============


Metadata
********

Author: **Unknown**


*Gridss was last updated on 2019-08-20*.
*This page was automatically generated on 2019-09-10*.
