# Parsing CWL

Common Workflow Language (CWL) is a common transport and execution format for workflows. Janis translates to CWL, and now Janis can *partially* translate from CWL.

## Quickstart

Make sure you have `janis-pipelines >= 0.11.1` AND `janis-pipelines.core >= 0.11.3` installed.

Note, we're using `janisdk` NOT the regular `janis` CLI:

```shell
janisdk fromcwl <yourworkflow.cwl>
```

If this translates correctly, this will produce python source code from your workflow.


## Case study

Result: https://github.com/PMCC-BioinformaticsCore/janis-pipelines/tree/master/janis_pipelines/kidsfirst

Converting Kids-First pipelines:

- [alignment workflow](https://github.com/kids-first/kf-alignment-workflow)
- [rnaseq workflow](https://github.com/kids-first/kf-rnaseq-workflow)
- [somatic workflow](https://github.com/kids-first/kf-somatic-workflow)
- [joint genotyping](https://github.com/kids-first/kf-jointgenotyping-workflow)

Clone one of the workflows to a location. Some things to watch out for:

- Ensure that no tool ids clash (eg: workflow.id, commandlinetool.id), Janis relies on these being unique IDS as they get used to generate variable names.
- Expressions that don't translate, these get logged as warnings.

For example, when translating the `kfdrc-rnaseq-workflow.cwl`, we find the following warnings:

```log
$ janisdk fromcwl /Users/franklinmichael/source/kf-rnaseq-workflow/workflow/kfdrc-rnaseq-workflow.cwl

[INFO]: Loading CWL file: /Users/franklinmichael/source/kf-rnaseq-workflow/workflow/kfdrc-rnaseq-workflow.cwl
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>"*TRIMMED." + inputs.readFilesIn1</expr>'
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>inputs.strand ? inputs.strand == "default" ? "" : "--"+inputs.strand : ""</expr>'
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>inputs.reads2 ? inputs.reads1.path+" "+inputs.reads2.path : "--single -l "+inputs.avg_frag_len+" -s "+inputs.std_dev+" "+inputs.reads1</expr>'
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>inputs.chimeric_sam_out.nameroot</expr>'
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>inputs.unsorted_bam.nameroot</expr>'
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>inputs.chimeric_sam_out.nameroot</expr>'
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>inputs.readFilesIn2 ? inputs.readFilesIn2.path : ''</expr>'
[WARN]: Couldn't translate javascript token, will use the placeholder '<expr>inputs.genomeDir.nameroot.split('.'</expr>'
[WARN]: Expression tools aren't well converted to Janis as they rely on unimplemented functionality: file:///Users/franklinmichael/source/kf-rnaseq-workflow/tools/expression_parse_strand_param.cwl#expression_strand_params```
```

## Limitations

- Can only convert very basic Javascript expressions 
- `when` conditionals aren't supported yet (we can't translate the expressions very well)
- Some requirements may not be faithfully represented
- Expression tools get translated, but they require extra Janis features to be implemented to work:
    - [Add ReadJsonOperator](https://github.com/PMCC-BioinformaticsCore/janis-core/issues/59)
    - [Allow Stdout / stderr to better participate in operations](https://github.com/PMCC-BioinformaticsCore/janis-core/issues/58)

Expressions are one of the trickiest bits of the translation, at the moment we can only parse some extremely basic expressions, for example:

- `${ return 800 }`: Literal value `800`
- `$(inputs.my_input)`: `InputSelector("my_input")`
- `$(inputs.my_input.basename)`: `BasenameOperator(InputSelector("my_input"))`
- `$(inputs.my_input.size)`: `FileSizeOperator(InputSelector("my_input"))`
- `$(inputs.my_input.contents)`: `ReadContents(InputSelector("my_input"))`


### Developer notes

Source code: [`janis_core/ingestion/fromcwl.py`](https://github.com/PMCC-BioinformaticsCore/janis-core/blob/master/janis_core/ingestion/fromcwl.py)

#### Expressions

Expressions are parsed using a combination of regular expressions and python string matching, there's nothing super clever there.

Tests: [`janis_core/tests/test_from_cwl.py`](https://github.com/PMCC-BioinformaticsCore/janis-core/blob/master/janis_core/tests/test_from_cwl.py)
