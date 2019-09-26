:orphan:


GATK4: MergeMutectStats
===============================================

Description
-------------

Tool identifier: ``GATK4MergeMutectStats``

Tool path: ``janis_bioinformatics.tools.gatk4.mergemutectstats.versions import Gatk4MergeMutectStats_4_1_2``

Version: 4.1.2.0

Container: ``broadinstitute/gatk:4.1.2.0``

Versions
*********

- `4.1.3.0 <gatk4mergemutectstats_4.1.3.0.html>`_
- 4.1.2.0 (current)

Documentation
-------------

URL
******
`TBD <TBD>`_

Tool documentation
******************
TBD

Outputs
-------
======  ========  ========================
name    type      documentation
======  ========  ========================
out     TextFile  Merged callability stats
======  ========  ========================

Inputs
------
Find the inputs below

Required inputs
***************

==========  ===============  ========  ==========  =================
name        type             prefix      position  documentation
==========  ===============  ========  ==========  =================
statsFiles  Array<TextFile>  --stats            0  Callability stats
==========  ===============  ========  ==========  =================

Optional inputs
***************

==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
mergedStatsOut  Optional<Filename>  -O                 1
==============  ==================  ========  ==========  ===============


Metadata
********

Author: Hollizeck Sebastian


*GATK4: MergeMutectStats was last updated on 2019-09-09*.
*This page was automatically generated on 2019-09-26*.
