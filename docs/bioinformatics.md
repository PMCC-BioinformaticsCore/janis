# Janis and Bioinformatics

Janis was designed with bioinformatics first because specifically cancer genomics was the scope it was built for. Janis is still applicable for lots of other domains, and shouldn't just be considered a bioinformatics workflow tool.

However there are a few special additions to Janis for bioinformatics research that we'd like to highlight.


## Indexed bams

Bam files often contain an index, often with the extension `.bai`, but the age old question is:

1. `myfile.bam.bai` (Janis preferred!)
2. `myfile.bai`

Most tools built with HTSLib accept both, but some pesky tools (looking at you [GATK](https://github.com/broadinstitute/gatk/issues/5299#issuecomment-565899113)) only accept `.bai`. 

Out of a small survey of our local bioinformaticians, we've chosen `.bam.bai` as the Janis preferred as we feel it creates a better link that the `.bai` index belongs to the bam file.

However, if your tool REQUIRES the use of this other format, you can do so with the following `secondaries_present_as` attribute:

The dictionary map is structured as: `{formatInType: desiredFormat}`

```python
# Tool input
ToolInput("myBam", BamBai, secondaries_present_as={".bai": "^.bai"})

# Tool output
ToolOutput("myOutput", BamBai, glob="file.bam", secondaries_present_as={".bai": "^.bai"})
```

