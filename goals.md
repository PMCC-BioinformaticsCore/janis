# Project Goals

## Introduction
There are more than 60 pipeline tools available and new tools appearing on an almost daily basis what is the
direction that the field is going in? The answer is that the complexity of the problem is also growing and
it becomes sensible to break the software problem into separate areas of concern. Historically, a single 
tool that could reliable run a few samples in a single address space may have been adequate. These days, we are
analysing large numbers of large samples with very complex workflows across very complex computational environments.

A recent initiative, Common Workflow Language (CWL, [commonwl.org](https://www.commonwl.org/)), is attempting to
create a standard for describing workflows. Other tools are starting to support this standard but there are also other
initiatives to separate the workflow definition from the workflow execution. The current state of play suggests the
following minimum decomposition of the problem space:

* Workflow definition languages such as CWL and Workflow Definition Language (WDL, 
[software.broadinstitute.org/wdl/](https://software.broadinstitute.org/wdl/)) from the Broad Institute
* Engines such as the Toil workflow engine (toil.readthedocs.io/en/latest/](https://toil.readthedocs.io/en/latest/) and 
the Broad Institute's Cromwell engine that will take a workflow definition
* Distributed resource managers such as legacy batch systems (Torque, slurm, etc) more modern solutions such as
mesos.


CWL, by necessity, is complex. It needs provide complete descriptions of:

* the job graph
* file handling
* software invocation

This detail is burdensome to provide. This exacerbated by the design of CWL which has limited programmatic abilities
for automating repetitive tasks. The design and documentation of CWL also makes working with it difficult. CWL is
probably more suitable as a transport format and publishable artefact than as primary user format for specifying
workflows.

There is also an opportunity for software solutions to help organisations provide organisational level support for their
users.

* To provide the ability to quickly develop pipelines for common use cases
* To provide appropriate clients for building, submitting, monitoring and managing workflows
* To devolve organisational best practice to uses. This is includes recommended
    * workflow steps
    * tools and tool settings
    * resource requirements
    * file handling


## Requirements
The requirements of this project are as follows:

Provide a mechanism for unsophisticated users to build text based workflow definitions. It is a requirement that it be
possible to specify the workflow in high level functional terms that have meaning to the data practitioners rather than
in terms of software tools, file dispositions, etc. To take a bioinformatics example, a simple variant calling pipeline
might look like:
```text
inputs:
  fastq:
    SequenceReadArchivePaired:
      forward-pattern: '*_R1.fastq.gz'
      backward-pattern: '*_R2.fastq.gz'
  ref:
    bam:
      path: 'path/to/reference'

steps:
  - step1:
      trim:
  - step2:
      align:
  - step3:
      dedup:
  - step4:
      call:
```

The solution can use a variety of mechanisms to complete the pipeline such as
* using textual order to imply execution order
* file types required for different steps can be used to determine dependencies
* convention

The implementation should emphasis ease of use over completeness. It is _not_ a requirement to be able to specify an
arbitrary job graph. Similarly, the graph will be known in advance, i.e. a dynamic job graph is not required. (CWL does 
not support it either.)

The implementation should provide a schema that can be consumed by other clients for rendering interfaces.

The implementation will provide an API so that tools and data types can be easily provided by the organisational unit
whose job it is to support pipelines and by advanced external uses. Other than examples and test cases, implementations
will not be part the project. The API will be as accessible as possible so, sadly, that means the API will be in Python.


## Michael's Cleanup Goals
- Add type annotations to everything
- Improve the logging format:
    - This will probably be a singleton class that can take logs from anywhere
    - Use levels for console and disk logs:
        - NONE (No log)
        - Critical (RED)
        - Info (White)
        - Debug (Grey)
    - Use color to distinguish levels
    - Write to file
    - Buffer Warnings and Errors
- Break the conversion into more discrete steps:
    - Load from YAML, JSON, XML, etc
    - Build DAG
- Ensure the graph nodes have the correct HASH function, and derive from some common (named) type
- Should labels across inputs / steps be unique? I feel like they should be unique....
- Remove the type inference from the workflow construction
    - Add warnings to user about matching types


 