
.. include:: performancevalidator_1.2.1

Performance Validator
============================================

Description
-------------

Tool identifier: ``performanceValidator``

Tool path: ``janis_bioinformatics.tools.validation.performancevalidator import PerformanceValidator_1_2_1``

Version: 1.2.1





Documentation
-------------

URL
******
*No URL to the documentation was provided*

Description
*********
*No documentation was provided: `contribute one <https://github.com/illusional>`_*

Outputs
-------
==================  ======  ===============
name                type    documentation
==================  ======  ===============
summaryMetrics      File
detailMetrics       File
contingencyMetrics  File
==================  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  ==========  ========  ==========  ===============
name       type        prefix    position    documentation
=========  ==========  ========  ==========  ===============
vcf        VCF
truth      VCFIDX
intervals  Array<VCF>
=========  ==========  ========  ==========  ===============

Optional inputs
***************

==========================  =================  ========  ==========  ===============
name                        type               prefix    position    documentation
==========================  =================  ========  ==========  ===============
treatMissingSitesAsHomeRef  Optional<Boolean>
==========================  =================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Performance Validator was last updated on **Unknown***.
*This page was automatically generated on 2019-07-24*.
