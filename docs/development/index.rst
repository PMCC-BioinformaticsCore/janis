Development
===========

|GitHub stars| |Build Status| |Documentation Status| |PyPI version|
|codecov|

    This project is work-in-progress and is still in developments.
    Although we welcome contributions, due to the immature state of this
    project we recommend raising issues through the `Github issues
    page <https://github.com/PMCC-BioinformaticsCore/janis/issues>`__.

Introduction
~~~~~~~~~~~~

This page has a collection of notes on the development of Janis.

Testing
-------

    *Further information*:
    `Testing <https://janis.readthedocs.io/en/latest/development/testing.html>`__

Within Janis there is a suite of code-level tests. High quality tests
allow fewer bugs to make it out to a production environment and gives
you a place to test your function.

All of the tests are placed in: ``janis/tests/test_*.py``. You can
create a file in that same directory with tests for your changes.

Running tests
~~~~~~~~~~~~~

In the terminal you can run the tests and generate a code coverage
report with the following bash command.

.. code:: bash

    nosetests -w janis --with-coverage --cover-package=janis

Releasing
---------

    *Further Information*:
    `Releasing <https://janis.readthedocs.io/en/latest/development/releasing.html>`__

Releasing is automatic! After you've run the tests, simply increment the
version number in ``setup.py`` (with resepect to
`SemVer <https://semver.org>`__), and tag that commit with the same
version identifier:

.. code:: bash

    git commit -m "Tag for v0.x.x release"
    git tag -a "v0.x.x" -m "Tag message"
    git push origin v0.x.x

Logger
------

This project uses a custom logger that can be imported from the root of
Janis:

.. code:: python

    from janis import Logger

-  ``Logger.mute()`` stops any output from reaching the console until
   ``unmute()`` is called. This does not affect logging to disk.
-  ``Logger.unmute()`` restores the logger's ability to write to the
   console if it has been muted. It will restore to the same *volume*
   before it was muted.
-  Logging functions:

   -  ``.log(str)`` - ``LogLevel.DEBUG``
   -  ``.info(str)`` - ``LogLevel.INFO``
   -  ``.warn(str)`` - ``LogLevel.WARNING``
   -  ``.critical(str)`` - ``LogLevel.CRITICAL``
   -  ``.log_exec(exception)`` - ``Loglevel.CRITICAL`` - Prepares
      exception in standard format

Adding a translation
--------------------

There are a few things that are super helpful to know before you add a
new translation. You should be familiar with Janis' representation of a
workflow.

    It's recommended you use an intermediary library (similar to cwlgen
    or python-wdlgen) to manage the representation. This allows you to
    focus on mapping concepts over syntax.

Concepts
~~~~~~~~

| There are a few things you'll need to keep in mind when you add a
translation, you should understand at least these workflow and janis
concepts:
|  - inputs - outputs - steps + tools - secondary files - command line
binding of tools - position, quoted - selectors - input selectors -
wildcard glob selectors - cpu and memory selectors - propagated resource
overrides

Translating
~~~~~~~~~~~

Add a file to the ``janis/translation`` with your new language's name.
Create a class called: ``$LangTranslator`` that should inherit from
``TranslatorBase``, and provide implementations for those methods.

Please keep your function sizes small and direct, and then write unit
tests to cover each component of the translation and then an integration
test of the whole translation on the related workflows.

You can find these in ``/janis/tests/test_translation_*.py``)

.. |GitHub stars| image:: https://img.shields.io/github/stars/PMCC-BioinformaticsCore/janis.svg?style=social
.. |Build Status| image:: https://travis-ci.org/PMCC-BioinformaticsCore/janis.svg?branch=master
   :target: https://travis-ci.org/PMCC-BioinformaticsCore/janis
.. |Documentation Status| image:: https://readthedocs.org/projects/janis/badge/?version=latest
   :target: https://janis.readthedocs.io/en/latest/?badge=latest
.. |PyPI version| image:: https://badge.fury.io/py/janis-pipelines.svg
   :target: https://badge.fury.io/py/janis-pipelines
.. |codecov| image:: https://codecov.io/gh/PMCC-BioinformaticsCore/janis/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/PMCC-BioinformaticsCore/janis