# Running the alignment

Running the workflow is the last piece of our tutorial, and we're going to use Janis to assist in this process.

## Requirements

In the [requirements](#), we specified that this tutorial required CWLTool and Docker to run the workflow.

## Using Janis-Runner

Janis will prepare a folder with your workflow and translated files at `~/janis/execution/$name/$date_$time_$tid/`. This is overridable with the `-o` parameter that must be an empty directory, Janis will create this output directory if it does not exist.

You can run the workflow with:

```bash
janis run alignment.py --engine cwltool
```

This will give you a status screen:
