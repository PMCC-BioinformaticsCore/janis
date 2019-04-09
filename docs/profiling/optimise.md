# Optimising bioinformatics tools

We have a number of of resource factors we can control to affect a tool's performance with regards to:
- Time taken
- Accuracy of results

For the most part, tool's scale with the amount of resources (more CPUs and more RAM), but so do the costs, so we want
to find a range in which these factors are balanced yet still give an appropriate running time.

In addition to the two metrics listed before, we include:
- Cost

as an additional factor in determining how we request resources for our tools, this is difficult to estimate
for local HPCs, however most [cloud](#cloud) vendors have fixed price lists based on the instance type and
the time the resource is used for.

### Reference files

For these tests, we use the _contigs renamed_ version of the human reference genome 
([HG38: NCBI](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/)).

The following table gives the file size information about the reference data sets we've used throughout this project,
including their secondary files using 
[CWL's secondary file syntax](https://www.commonwl.org/v1.0/CommandLineTool.html#CommandInputParameter).

| Name                  | Size    | Secondary Files                            |
|-----------------------|---------|--------------------------------------------|
|              Assembly | 5.34 GB | .fai, .amb, .ann, .bwt,  .pac, .sa, ^.dict |
|          SNPS (DBSNP) | 1.56 GB | .tbi                                       |
|         SNPS (1000GP) | 1.89 GB | .tbi                                       |
|          Known Indels | 63.3 MB | .tbi                                       |
| Mills Indels (1000Gp) | 22.1 MB | .tbi                                       |


## Cloud support

Talk about:
- Optimising for specific instance types (grouped resources can save money over custom types)


## Other technologies:

- FPGAs (Field programmable gate arrays)
- TPUs (Tensor processing units)


## How Janis works

For some tools, we have gone through the process of determining a happy medium between time taken to run the task
for each capture type vs the cost of running on cloud vendors. For the most part, we've tried to keep the values on
cloud and local similar, but as discussed in the [cloud](#cloud) section, there may be some discrepancies.

These values are provided within the tool definition, and a resource file can be produced for specific vendors that
may be used when running to try and optimise the workflow running experience.

```python
w = MyWorkflow()      # My constructed workflow
w.translate("cwl", write_resource_file="gcp")
```
