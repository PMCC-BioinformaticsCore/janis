:orphan:


Concatenate
=================

Description
-------------

Tool identifier: ``cat``

Tool path: ``janis_unix.tools.cat import Cat``

Version: v1.0.0

Container: ``ubuntu:latest``



Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
The cat utility reads files sequentially, writing them to the standard output. The file operands are processed in command-line order. If file is a single dash (`-') or absent,cat reads from the standard input. If file is a UNIX domain socket, cat connects to it and then reads it until EOF. This complements the UNIX domain binding capability available in inetd(8).

Outputs
-------
======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>
======  ============  ===============

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


*Concatenate was last updated on 2019-07-26 00:00:00*.
*This page was automatically generated on 2019-09-10*.
