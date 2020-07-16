:orphan:

VcfLib: Vcf Filter
==============================

*1 contributor · 1 version*

Filter the specified vcf file using the set of filters.
Filters are specified in the form "<ID> <operator> <value>:
 -f "DP > 10"  # for info fields
 -g "GT = 1|1" # for genotype fields
 -f "CpG"  # for 'flag' fields

Operators can be any of: =, !, <, >, |, &

Any number of filters may be specified.  They are combined via logical AND
unless --or is specified on the command line.  Obtain logical negation through
the use of parentheses, e.g. "! ( DP = 10 )"

For convenience, you can specify "QUAL" to refer to the quality of the site, even
though it does not appear in the INFO fields.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcffilter.versions import VcfFilter_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcffilter_step",
           VcfFilter_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcffilter_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcffilter:

.. code-block:: bash

   # user inputs
   janis inputs vcffilter > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run vcffilter with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcffilter





Information
------------


:ID: ``vcffilter``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Michael Franklin
:Citations: None
:Created: 2020-06-04
:Updated: 2020-06-04



Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<VCF>  Filtered VCF
======  ===========  ===============



Additional configuration (inputs)
---------------------------------

===============  =========================  =================  ==========  ===================================================================================================================================================================
name             type                       prefix               position  documentation
===============  =========================  =================  ==========  ===================================================================================================================================================================
vcf              VCF                                                    1  VCF to filter
info_filter      Optional<String>           --info-filter                  (-f) specifies a filter to apply to the info fields of records, removes alleles which do not pass the filter
genotype_filter  Optional<String>           --genotype-filter              (-g) specifies a filter to apply to the genotype fields of records
keep_info        Optional<Boolean>          --keep-info                    (-k) used in conjunction with '-g', keeps variant info, but removes genotype
filter_sites     Optional<Boolean>          --filter-sites                 (-s) filter entire records, not just alleles
tag_pass         Optional<String>           --tag-pass                     (-t) tag vcf records as positively filtered with this tag, print all records
tag_fail         Optional<String>           --tag-fail                     (-F) tag vcf records as negatively filtered with this tag, print all records
append_filter    Optional<Boolean>          --append-filter                (-A) append the existing filter tag, don't just replace it
allele_tag       Optional<String>           --allele-tag                   (-a) apply -t on a per-allele basis. adds or sets the corresponding INFO field tag
invert           Optional<Boolean>          --invert                       (-v) inverts the filter, e.g. grep -v
use_logical_or   Optional<Boolean>          --or                           (-o) use logical OR instead of AND to combine filters
region           Optional<Array<BedTABIX>>  --region                       (-r) specify a region on which to target the filtering, requires a BGZF compressed file which has been indexed with tabix.  any number of regions may be specified.
===============  =========================  =================  ==========  ===================================================================================================================================================================
