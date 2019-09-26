:orphan:


GATK4: LearnReadOrientationModel
=================================================================

Description
-------------

Tool identifier: ``GATK4LearnReadOrientationModel``

Tool path: ``janis_bioinformatics.tools.gatk4.learnreadorientationmodel.versions import Gatk4LearnReadOrientationModel_4_1_2``

Version: 4.1.2.0

Container: ``broadinstitute/gatk:4.1.2.0``

Versions
*********

- `4.1.3.0 <gatk4learnreadorientationmodel_4.1.3.0.html>`_
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
======  =================  =============================================================
name    type               documentation
======  =================  =============================================================
out     CompressedTarFile  Model file containing information about fragment orientations
======  =================  =============================================================

Inputs
------
Find the inputs below

Required inputs
***************

===============  ========================  ========  ==========  ============================================
name             type                      prefix      position  documentation
===============  ========================  ========  ==========  ============================================
f1r2CountsFiles  Array<CompressedTarFile>  -I                 0  Counts for the read orientation of fragments
===============  ========================  ========  ==========  ============================================

Optional inputs
***************

===============  ==================  ===================  ==========  =======================================================
name             type                prefix                 position  documentation
===============  ==================  ===================  ==========  =======================================================
numEmIterations  Optional<Integer>   --num-em-iterations           1  Amount of iterations for the em process before it bails
modelFileOut     Optional<Filename>  -O                            3
===============  ==================  ===================  ==========  =======================================================


Metadata
********

Author: Hollizeck Sebastian


*GATK4: LearnReadOrientationModel was last updated on 2019-09-09*.
*This page was automatically generated on 2019-09-26*.
