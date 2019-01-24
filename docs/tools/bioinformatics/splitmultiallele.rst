
splitmultiallele
================
*bioinformatics*

Documentation
-------------

Docker
******
``heuermh/vt``

URL
******
*No URL to the documentation was provided*: `contribute one <https://github.com/illusional>`_

Docstring
*********
VcfSplitMultiAllele.sh
    
    Currently stored at: '/researchers/jiaan.yu/WGS_pipeline/VcfSplitMultiAllele.sh’
    
    This CommandTool is an attempt to translate the shell script as a CommandTool.
    
    It uses the commands 'sed' and 'vt' (where $0: input, $1: output) as the following:
    
        > sed 's/ID=AD,Number=./ID=AD,Number=R/' < $1       |            vt decompose -s - -o -                          |            vt normalize -n -q - -o - -r $HumanREF          |            sed 's/ID=AD,Number=./ID=AD,Number=1/' > $2
        
    SED documentation:
        Usage: sed [OPTION]... {script-only-if-no-other-script} [input-file]...
        
          -n, --quiet, --silent
                         suppress automatic printing of pattern space
          -e script, --expression=script
                         add the script to the commands to be executed
          -f script-file, --file=script-file
                         add the contents of script-file to the commands to be executed
          --follow-symlinks
                         follow symlinks when processing in place
          -i[SUFFIX], --in-place[=SUFFIX]
                         edit files in place (makes backup if SUFFIX supplied)
          -c, --copy
                         use copy instead of rename when shuffling files in -i mode
          -b, --binary
                         does nothing; for compatibility with WIN32/CYGWIN/MSDOS/EMX (
                         open files in binary mode (CR+LFs are not treated specially))
          -l N, --line-length=N
                         specify the desired line-wrap length for the `l' command
          --posix
                         disable all GNU extensions.
          -r, --regexp-extended
                         use extended regular expressions in the script.
          -s, --separate
                         consider files as separate rather than as a single continuous
                         long stream.
          -u, --unbuffered
                         load minimal amounts of data from the input files and flush
                         the output buffers more often
          -z, --null-data
                         separate lines by NUL characters
          --help
                         display this help and exit
          --version
                         output version information and exit
        
        If no -e, --expression, -f, or --file option is given, then the first
        non-option argument is taken as the sed script to interpret.  All
        remaining arguments are names of input files; if no input files are
        specified, then the standard input is read.
        
        GNU sed home page: <http://www.gnu.org/software/sed/>.
        General help using GNU software: <http://www.gnu.org/gethelp/>.
        
    VT decompose documentation:
        options : -s  smart decomposition [false]
              -d  debug [false]
              -f  filter expression []
              -o  output VCF file [-]
              -I  file containing list of intervals []
              -i  intervals []
              -?  displays help
              
    VT normalize documentation:
        options : -o  output VCF file [-]
              -d  debug [false]
              -q  do not print options and summary [false]
              -m  warns but does not exit when REF is inconsistent
                  with masked reference sequence for non SNPs.
                  This overides the -n option [false]
              -n  warns but does not exit when REF is inconsistent
                  with reference sequence for non SNPs [false]
              -f  filter expression []
              -w  window size for local sorting of variants [10000]
              -I  file containing list of intervals []
              -i  intervals []
              -r  reference sequence fasta file []
              -?  displays help

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
output  VCF
======  ======  ===============

Inputs
------
==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
input           VCF                                    2
reference       FastaWithDict       -r                 7
outputFilename  Optional<Filename>  >                 10
==============  ==================  ========  ==========  ===============


*This page was automatically generated*
