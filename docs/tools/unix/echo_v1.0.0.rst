:orphan:


Echo
===========

Description
-------------

Tool identifier: ``echo``

Tool path: ``janis_unix.tools.echo import Echo``

Version: v1.0.0

Container: ``ubuntu:latest``



Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
The echo utility writes any specified operands, separated by single blank (` ') characters and followed by a newline (`
') character, to the standard output.

Some shells may provide a builtin echo command which is similar or identical to this utility. Most notably, the builtin echo in sh(1) does not accept the -n option. Consult the builtin(1) manual page.

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
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
inp     String                     1
======  ======  ========  ==========  ===============

Optional inputs
***************

===============  =================  ========  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================
name             type               prefix    position    documentation
===============  =================  ========  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================
include_newline  Optional<Boolean>  -n                    Do not print the trailing newline character.  This may also be achieved by appending `\c' to the end of the string, as is done by iBCS2 compatible systems.  Note that this option as well as the effect of `\c' are implementation-defined in IEEE Std 1003.1-2001 (``POSIX.1'') as amended by Cor. 1-2002.  Applications aiming for maximum portability are strongly encouraged to use printf(1) to suppress the newline character.
===============  =================  ========  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: **Unknown**


*Echo was last updated on **Unknown***.
*This page was automatically generated on 2019-09-26*.
