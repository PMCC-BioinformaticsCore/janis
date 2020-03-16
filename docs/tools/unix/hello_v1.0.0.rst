:orphan:

Hello, World!
=====================

*1 contributor · 1 version*

:ID: ``hello``
:Python: ``janis_unix.tools.hello import HelloWorkflow``
:Versions: v1.0.0
:Authors: Michael Franklin
:Citations: 
:Created: None
:Updated: 2019-08-12
:Required inputs:

:Outputs: 
   - ``out: File``

Documentation
-------------

URL: *No URL to the documentation was provided*

This is the 'Hello, world' equivalent workflow that uses the Echo unix
tool to log "Hello, World!" to the console, and collects the result.

This is designed to be the first example that you can run with janis, ie:
    
``janis run hello``


Embedded Tools
***************

====  ===============
Echo  ``echo/v1.0.0``
====  ===============

------

Additional configuration (inputs)
---------------------------------

======  ================  ===============
name    type              documentation
======  ================  ===============
inp     Optional<String>
======  ================  ===============

.
