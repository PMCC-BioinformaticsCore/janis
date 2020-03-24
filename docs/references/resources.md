# Configuring resources (CPU / Memory)

Sometimes you'll want to override the amount of resources a tool will get.


## Hints

> There are changes to the tool docs coming to encapsulate this information.

Some tools are aware of hints, and can change their resources based on the hints you provide. For example, BWA MEM responds to the `captureType` hint (`targeted`, `exome`, `chromosome`, `30x`, `90x`).


## Generating resources template

Sometimes you want complete custom control over which resources each of your tool, you can generate an inputs template for resources with the following:

```bash
$ janis inputs --resources BWAAligner

# bwamem_runtime_cpu: 16
# bwamem_runtime_disks: local-disk 60 SSD
# bwamem_runtime_memory: 16
# cutadapt_runtime_cpu: 5
# cutadapt_runtime_disks: local-disk 60 SSD
# cutadapt_runtime_memory: 4
# fastq: null
# reference: null
# sample_name: null
# sortsam_runtime_cpu: 1
# sortsam_runtime_disks: local-disk 60 SSD
# sortsam_runtime_memory: 8
```

## Limitations

- You are unable to size a specific shard within a scatter. (See _Expressions in runtime attributes_) 
- The only way to set the `disk` size for a tool is to generate the inputs template, and set a value.
- Currently (2020-03-20) you are unable to change the time limit of tasks

## Upcoming work

> Leave an [issue on GitHub](https://github.com/PMCC-BioinformaticsCore/janis-assistant/issues/new) (or comment on the linked issue) if you have comments.

_Proposed_: [PMCC-BioinformaticsCore/janis-assistant#9](https://github.com/PMCC-BioinformaticsCore/janis-assistant/issues/9)


- Providing hints to get appropriately sized resources

    ```bash
    janis inputs [PROPOSED: --hint-captureType 30x] --resource WGSGermlineGATK
    ```

_Proposed_: 

- Expressions in runtime attributes.