Unit Test Framework
========================

.. note::
	Available in v0.11.0 and later

Overview
*********

You can write test cases for your tool by defining the ``tests()`` function in your Tool class.
Test cases defined by this function will be picked up by ``janisdk run-test`` command.

This unit test framework provides several predefined preprocessors to
transform Janis execution output data in the format that can be tested.
For example, there is a preprocessor to read the md5 checksum of an output file.
In addition to the predefined preprocessors, the frameworks also allows users to define and pass their own preprocessors.


Define test cases
-----------------
You can define multiple test cases per tool.
For each test case, you can declare multiple expected outputs.
Mostly, you really only need to define more test cases if they require different input data.

The classes we use here :class:`janis_core.tool.test_classes.TTestCase`, :class:`janis_core.tool.test_classes.TTestExpectedOutput` and :class:`janis_core.tool.test_classes.TTestPreprocessor` are declared at the bottom of this document.

.. code-block:: python

    class BwaAligner(BioinformaticsWorkflow):
        def id(self):
            return "BwaAligner"

        ...

        def tests(self):
            return [
                TTestCase(
                    name="basic",
                    input={
                        "bam": "https://some-public-container/directory/small.bam"
                    },
                    output=[
                        TTestExpectedOutput(
                            tag="out",
                            preprocessor=TTestPreprocessor.FileMd5,
                            operator=operator.eq,
                            expected_value="dc58fe92a9bb0c897c85804758dfadbf",
                        ),
                        TTestExpectedOutput(
                            tag="out",
                            preprocessor=TTestPreprocessor.FileContent,
                            operator=operator.contains,
                            expected_value="19384 + 0 in total (QC-passed reads + QC-failed reads)",
                        ),
                        TTestExpectedOutput(
                            tag="out",
                            preprocessor=TTestPreprocessor.LineCount,
                            operator=operator.eq,
                            expected_value=13,
                        ),
                    ],
                )
            ]

Run the tests
-------------
.. code-block:: console

    # Run a specific test case
    janisdk run-test --test-case [TEST CASE NAME] [TOOL ID]

    # Run ALL test cases of one tool
    janisdk run-test [TOOL ID]

To run the example test case shown above:

.. code-block:: console

    # Run a specific test case
    janisdk run-test --test-case=basic BwaAligner

    # Run all test cases
    janisdk run-test BwaAligner

Test Files
**********

There are two different ways to store your test files (input and expected output files):

Remote HTTP files:
------------------

you can use a publicly accessible http link ``https://some-public-container/directory/small.bam``.

* Input files will be downloaded to a cache folder in ``~/.janis/remote_file_cache`` folder. This is the same directory where files will be cached when you run ``janis run``.
* Expected output files however will be cached in the test directory ``[WORKING DIRECTORY WHERE janisdk run-test is run]/tests_output/cached_test_files/``.

If the same url is found in the cache directory, we will not re-download the files unless the ``Last-Modified`` http header has changed. If you want to force the files to be re-downloaded, you will need to remove the files from the cache directories.

Local test files:
----------------

you can store your files in local directory named ``test_data``. There are a few different examples of where you can place this directory.
Example from ``janis-bioinformatics`` project:

* A ``test_data`` folder that contain files to be shared by multiple tools can be located at ``janis_bioinformatics/tools/test_data``. To access files in this directory, you can call ``os.path.join(BioinformaticsTool.test_data_path(), "small.bam")``.
* A ``test_data`` folder that contain files to be used by ``flagstat`` can be located at ``janis_bioinformatics/tools/samtools/flagstat/test_data``. To access files in this directory from within the ``SamToolsFlagstatBase`` class, you can call ``os.path.join(self.test_data_path(), "small.bam")``.

Preprocessors and Comparison Operators
**************************************

``TTestExpectedOutput.preprocessor`` is used to reformat the Tool output.
``TTestExpectedOutput.operator`` is used to compare output value with the expected output.

Predefined Preprocessors
------------------------

* **Value**: No preprocessing, value as output by Janis e.g. an integer, string, or a file path for a File type output.
* **FileDiff**: The differences between two files as output by ``difflib.unified_diff``. This can only be applied on File type output. If this preprocessor is used, ``TTestExpectedOutput.file_diff_source`` must be provided. ``file_diff_source`` must contain the file path to compare the output file with.
* **LinesDiff** Number of different lines between two files as a tuple ``(additions, deletions)``, as diff'd by the ``FileDiff`` preprocessor. If this preprocessor is used, ``TTestExpectedOutput.file_diff_source`` must be provided. ``file_diff_source`` must contain the file path to compare the output file with.
* **FileContent**: Extract the file content. This can only be applied to File type output.
* **FileExists**: Check if a file exists. It returns a True/False value. This can only be applied to File type output.
* **FileSize**: File size is bytes. This can only be applied on File type output.
* **FileMd5**: Md5 checksum of a file. This can only be applied to File type output.
* **LineCount**: Count the number of lines in a string or in a file.
* **ListSize**: Count the number of items in a list. This can only be applied to Array type output.


Custom preprocessor example:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example below, we are testing a tool that has an output field named ``out``.
The output value of this field is a file path that points to the location of a BAM file.
We want to test the flagstat value of this BAM file.
Here, we define your custom preprocessor function that
takes a file path as input and returns a string that contains the flagstat value of a BAM file.


``TTestExpectedOutput.expected_file`` simply points to a file that contains the expected output value.
You can also replace this with ``TTestExpectedOutput.expected_value="19384 + 0 in total (QC-passed reads + QC-failed reads)\n ..."``

.. code-block:: python

    TTestExpectedOutput(
        tag="out",
        preprocessor=Bam.flagstat,
        operator=operator.eq,
        expected_file="https://some-public-container/directory/flagstat.txt"
    )

.. code-block:: python

    class Bam(File):
        ...

        @classmethod
        def flagstat(cls, file_path: str):
            command = ["samtools", "flagstat", file_path]
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )

            if result.stderr:
                raise Exception(result.stderr)

            return result.stdout


Custom operator example:
^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we also want to test the flagstat value of BAM file returned by the ``out`` output field.

Here, instead of writing a custom preprocessors, we write a custom operator that takes two file path and compare the flagstat output of this two files.

.. code-block:: python

    TTestExpectedOutput(
        tag="out",
        preprocessor=TTestPreprocessor.Value,
        operator=Bam.equal,
        expected_value="https://some-public-container/directory/small.bam"
    )

.. code-block:: python

    class Bam(File):
        ...

        @classmethod
        def equal(cls, file_path_1: str, file_path_2: str):
            flagstat1 = cls.flagstat(file_path_1)
            flagstat2 = cls.flagstat(file_path_2)

            return flagstat1 == flagstat2

Declaration
***********

.. autoclass:: janis_core.tool.test_classes.TTestCase

.. autoclass:: janis_core.tool.test_classes.TTestExpectedOutput

.. autoclass:: janis_core.tool.test_classes.TTestPreprocessor


