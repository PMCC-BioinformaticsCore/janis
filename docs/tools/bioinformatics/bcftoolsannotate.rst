
bcftoolsannotate
================
*bioinformatics*

Documentation
-------------

URL
******
`https://samtools.github.io/bcftools/bcftools.html#annotate <https://samtools.github.io/bcftools/bcftools.html#annotate/>`_

Docstring
*********
BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary 
    counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed.

    Most commands accept VCF, bgzipped VCF and BCF with filetype detected automatically even when streaming 
    from a pipe. \Indexed VCF and BCF will work in all situations. Un-indexed VCF and BCF and streams will 
    work in most, but not all situations. In general, whenever multiple VCFs are read simultaneously, 
    they must be indexed and therefore also compressed.
    
    BCFtools is designed to work on a stream. It regards an input file "-" as the standard input (stdin) 
    and outputs to the standard output (stdout). Several commands can thus be combined with Unix pipes.

    Documentation: https://samtools.github.io/bcftools/bcftools.html
    -------------------------------------------------------------------------

    Add or remove annotations.
    

*This page was automatically generated*
