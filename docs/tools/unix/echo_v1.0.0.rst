:orphan:

Echo
===========

0 contributors · 1 version

:ID: ``echo``
:Python: ``janis_unix.tools.echo import Echo``
:Versions: v1.0.0
:Container: ubuntu:latest
:Authors: 
:Citations: None
:Created: None
:Updated: None
:Required inputs:
   - ``inp: String``
:Outputs: 
   - ``out: stdout<File>``

Documentation
-------------

URL: *No URL to the documentation was provided*

The echo utility writes any specified operands, separated by single blank (` ') characters and followed by a newline (`
') character, to the standard output.

Some shells may provide a builtin echo command which is similar or identical to this utility. Most notably, the builtin echo in sh(1) does not accept the -n option. Consult the builtin(1) manual page.

------

Additional configuration (inputs)
---------------------------------

===============  =================  =====================================================================================================================================================================================================================================================================================================================================================================================================================================
name             type               documentation
===============  =================  =====================================================================================================================================================================================================================================================================================================================================================================================================================================
inp              String
include_newline  Optional<Boolean>  Do not print the trailing newline character.  This may also be achieved by appending `\c' to the end of the string, as is done by iBCS2 compatible systems.  Note that this option as well as the effect of `\c' are implementation-defined in IEEE Std 1003.1-2001 (``POSIX.1'') as amended by Cor. 1-2002.  Applications aiming for maximum portability are strongly encouraged to use printf(1) to suppress the newline character.
===============  =================  =====================================================================================================================================================================================================================================================================================================================================================================================================================================

