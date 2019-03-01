Workflow Specifications
***********************

*This page is under construction*

Workflows can be specified using many different languages, in fact there are `almost 150 pipeline toolkits <https://github.com/pditommaso/awesome-pipeline>`_.

The main two specifications we focus on here are:

1. `Common Workflow Language <https://github.com/common-workflow-language/common-workflow-language>`_ (CWL)
2. `Workflow Description Language <https://software.broadinstitute.org/wdl/>`_ (WDL)

One of the main goals of Janice is to avoid writing a workflow in a single specification, but rather in a generic format which can be `transpiled <https://www.stevefenton.co.uk/2012/11/compiling-vs-transpiling/>`_ into CWL or WDL. We also get the advantage of adding custom types that make connecting pipelines together and allow you to forget about pesky secondary files.


How to specify workflows in Janis
=================================

There are 2 methods to specifying a Janis worfklow:

1. Use the Python :mod:`Pipeline` classes, see a `user guide here </userguide.html>`_.
2. Use the WEHI simplified workflow specification language.



