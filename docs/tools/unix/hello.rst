
.. include:: hello_v0.1.0

Hello, World!
=====================

Description
-------------

Tool identifier: ``hello``

Tool path: ``janis_unix.tools.hello import HelloWorkflow``

Version: v0.1.0





Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
This is the 'Hello, world' equivalent workflow that uses the Echo unix
tool to log "Hello, World!" to the console, and collects the result.

This is designed to be the first example that you can run with janis, ie:
    
``janis run hello``


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

======  ================  ========  ==========  ===============
name    type              prefix    position    documentation
======  ================  ========  ==========  ===============
inp     Optional<String>
======  ================  ========  ==========  ===============


Metadata
********

Author: **Unknown**


*Hello, World! was last updated on 2019-08-12*.
*This page was automatically generated on 2019-09-10*.
