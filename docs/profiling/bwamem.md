# Profiling BWA Mem

BWA Mem is a common aligner that produces a SAM file, from Fastq pair.

We're using the following version of BWA Mem:

- ``

Our inputs are the human reference genome and the FastQ pair produced 
by `cutadapt` (along with the headerGroupString (`-r`)). 
The threads option (`-t`) is scaled 

## Results

| CPU \ Mem |  8GB  | 12GB | 14GB  | 
|-----------|-------|------|-------| 
| 4 cores   | 1483s |  -   | 1594s | 
| 8 cores   | 884s  |  -   | 846s  | 
| 12 cores  | 661s  |  -   | 659s  | 
| 16 cores  | 434s  |  -   |   -   | 
| 20 cores  | 350s  | 350s |   -   | 
| 24 cores  | _DNF_ | 256s |   -   | 
| 28 cores  | _DNF_ | 236s |   -   | 

Interestingly, BWA mem doesn't perform significantly better when there is more ram given.
The exception to this is as the cores scale, the task might not complete without a
sufficient amount of RAM. In the diagram below, we plot the number of cores vs time to complete:



In this example, we have a lot more missing data points, primarily due to the longer
running task (between 6-30 minutes) for each run.

Similar to `cutadapt`, we can compute the approximate CPU time per job, 
and see how BWA scales with number of processes.

|  Cores   | 8GB  | 12-14GB |
|----------|------|---------|
| 4 cores  | 5932 | 6376    |
| 8 cores  | 7072 | 6768    |
| 12 cores | 7932 | 7908    |
| 16 cores | 6944 |   -     |
| 20 cores | 7000 | 7000    |
| 24 cores |_DNF_ | 6144    |
| 28 cores |_DNF_ | 6608    |

In some cases, the CPU improves as it is parallelised. This might be a side effect of the 
batch scheduler giving the job different hardware as the number of required processors grew.

## Cost

> Although it's not strictly correct as the platform might decide to choose a more efficient analysis and
costs scale in different ways, we can get a decent approximation of the compute cost 
by using the following custom machine values (for early 2019, Australia):
>
> | Type | [GCP](https://cloud.google.com/compute/pricing#custommachinetypepricing) |
> |------|---------------------------|
> | CPU  |  $ 1.24667E-05 / (cpu-s)  |
> | RAM  | $ 1.66944E-06 / (gig-sec) |

We can generate the following cost table: 

| Cost | 8          | 12         | 14         | 
|------|------------|------------|------------| 
| 4    |  $0.09376  |            |  $0.11674  | 
| 8    |  $0.09997  |            |  $0.10415  | 
| 12   |  $0.10771  |            |  $0.11399  | 
| 16   |  $0.09236  |            | -          | 
| 20   |  $0.09194  |  $0.08727  |            | 
| 24   |            |  $0.07660  |            | 
| 28   |            |  $0.08238  |            | 

When we plot this, we see a slightly irregular cost graph, especially given how the number of cores
vs total time taken to process approximately trended `log2`.

Our costs are approximately minimized on a chromosomal level at approximately 24 cores.
BWA mem doesn't perform substantially better with more RAM, and the minimum memory required to
run with these number of cores was 12GB.

