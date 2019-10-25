:orphan:

Performance Validator
============================================

0 contributors · 1 version

:ID: ``performanceValidator``
:Python: ``janis_bioinformatics.tools.validation.performancevalidator import PerformanceValidator_1_2_1``
:Versions: 1.2.1
:Authors: 
:Citations: 
:Created: None
:Updated: None
:Required inputs:
   - ``vcf: VCF``

   - ``truth: IndexedVCF``

   - ``intervals: Array<VCF>``
:Outputs: 
   - ``summaryMetrics: File``

   - ``detailMetrics: File``

   - ``contingencyMetrics: File``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

===========================  ====================================
BGZip                        ``bgzip/1.2.1``
Tabix                        ``tabix/1.2.1``
GATK4: Genotype Concordance  ``Gatk4GenotypeConcordance/4.1.4.0``
===========================  ====================================

------

Additional configuration (inputs)
---------------------------------

==========================================  =================  ===============
name                                        type               documentation
==========================================  =================  ===============
vcf                                         VCF
truth                                       IndexedVCF
intervals                                   Array<VCF>
genotypeConcord_treatMissingSitesAsHomeRef  Optional<Boolean>
==========================================  =================  ===============

.
