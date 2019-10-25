:orphan:

GATK4: LearnReadOrientationModel
=================================================================

1 contributor · 3 versions

:ID: ``GATK4LearnReadOrientationModel``
:Python: ``janis_bioinformatics.tools.gatk4.learnreadorientationmodel.versions import Gatk4LearnReadOrientationModel_4_1_3``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09
:Required inputs:
   - ``f1r2CountsFiles: Array<CompressedTarFile>``
:Outputs: 
   - ``out: CompressedTarFile``

Documentation
-------------

URL: `TBD <TBD>`_

TBD

------

Additional configuration (inputs)
---------------------------------

===============  ========================  =======================================================
name             type                      documentation
===============  ========================  =======================================================
f1r2CountsFiles  Array<CompressedTarFile>  Counts for the read orientation of fragments
numEmIterations  Optional<Integer>         Amount of iterations for the em process before it bails
modelFileOut     Optional<Filename>
===============  ========================  =======================================================

