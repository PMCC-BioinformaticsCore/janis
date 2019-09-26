:orphan:


Merge and Mark Duplicates
============================================

Description
-------------

Tool identifier: ``mergeAndMarkBams``

Tool path: ``janis_bioinformatics.tools.common.mergeandmark.mergeandmark_4_0 import MergeAndMarkBams_4_0``

Version: 4.0.12



Versions
*********

- `4.1.3 <mergeandmarkbams_4.1.3.html>`_
- 4.0.12 (current)

Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
======  =======  ===============
name    type     documentation
======  =======  ===============
out     BamPair
======  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ==============  ========  ==========  ===============
name    type            prefix    position    documentation
======  ==============  ========  ==========  ===============
bams    Array<BamPair>
======  ==============  ========  ==========  ===============

Optional inputs
***************

==================================  =================  ========  ==========  ===============
name                                type               prefix    position    documentation
==================================  =================  ========  ==========  ===============
createIndex                         Optional<Boolean>
maxRecordsInRam                     Optional<Integer>
mergeSamFiles_useThreading          Optional<Boolean>
mergeSamFiles_validationStringency  Optional<String>
==================================  =================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Merge and Mark Duplicates was last updated on **Unknown***.
*This page was automatically generated on 2019-09-26*.
