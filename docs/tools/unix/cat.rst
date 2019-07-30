
.. include:: cat_latest

Concatenate
=================

Description
-------------

Tool identifier: ``cat``

Tool path: ``janis_unix.tools.cat import Cat``

Version: latest

Container: ``ubuntu:latest``



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
======  ======  ===============
name    type    documentation
======  ======  ===============
out     Stdout
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ======  ========  ==========  ===============
name    type    prefix    position    documentation
======  ======  ========  ==========  ===============
======  ======  ========  ==========  ===============

Optional inputs
***************

====================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================
name                  type                   prefix      position  documentation
====================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================
file                  Optional<File>
files                 Optional<Array<File>>                     1
numberOutput          Optional<Boolean>      -n                    Number the output lines, starting at 1.
numberNonBlank        Optional<Boolean>      -b                    Number the non-blank output lines, starting at 1.
disableOutputBuffer   Optional<Boolean>      -u                    Disable output buffering.
squeeze               Optional<Boolean>      -s                    Squeeze multiple adjacent empty lines, causing the output to be single spaced.
displayNonprintChars  Optional<Boolean>      -e                    Display non-printing characters (see the -v option), and display a dollar sign (`$') at the end of each line.
displayNon            Optional<Boolean>      -t                    Display non-printing characters (see the -v option), and display tab characters as `^I'.
numberNonBlank        Optional<Boolean>      -v                    Display non-printing characters so they are visible.  Control characters print as `^X' for control-X; the delete character (octal 0177) prints as `^?'.  Non-ASCII characters (with the high bit set) are printed as `M-' (for meta) followed by the character for the low 7 bits.
====================  =====================  ========  ==========  ==================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*Concatenate was last updated on **Unknown***.
*This page was automatically generated on 2019-07-30*.
