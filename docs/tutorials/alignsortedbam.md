
# Building a simple bioinformatics workflow  
  
Welcome to this tutorial on building a simple bioinformatics workflow! Prior to starting this, we recommend completing these tutorials:  
  
- [Getting started with Janis](https://janis.readthedocs.io/en/latest/tutorials/gettingstarted.html)  
- [Constructing your first workflow](https://janis.readthedocs.io/en/latest/tutorials/simple.html)  
  
## Overview  
  
In this tutorial we're going to build a simple bioinformatics workflow that converts a `Fastq` (sequence data), into an aligned `Bam + Bai` pair (aligned sequence data). This is the [`AlignSortedBam`](https://janis.readthedocs.io/en/tools/bioinformatics/alignsortedbam.html) workflow packaged with the bioinformatics tool shed.  
  
To build our simple bioinformatics workflow, we're going to:  
  
- Subclass the `janis.Workflow` class,  
- Import the Tools and DataTypes from the Janis Bioinformatics toolshed,  
- Declare a `janis.Step` for each step in our pipeline,  
- Connect inputs and steps together  
- Get our output  
  
### What you'll need  
  
> Further info: [Getting Started](https://janis.readthedocs.io/en/latest/tutorials/gettingstarted.html)  
  
To complete this tutorial, you'll need to have an full installation of janis (see the [Getting Started](https://janis.readthedocs.io/en/latest/tutorials/gettingstarted.html) guide for more info). You can install janis with the bioinformatics tools by running:  
```bash  
pip install janis-pipelines[bioinformatics]  
```  
  
## Picking the tools  
  
For this workflow, we're going to create a 4-step Fastq alignment pipeline, using the following tools (in this order):  
  
1. [Bwa Mem](https://janis.readthedocs.io/en/latest/tools/bioinformatics/bwa/bwamem.html#inputs)  
2. [Samtools View](https://janis.readthedocs.io/en/latest/tools/bioinformatics/samtools/samtoolssort.html#inputs)  
3. [GATK4 SortSam](https://janis.readthedocs.io/en/latest/tools/bioinformatics/GATK4/gatk4sortsam.html#inputs)  
  
These three tools will run in sequence (the next tools depends on the previous), and there is no scatter / gathering involved in this pipeline to keep things simple.  
  
It's recommended you visit the documentation pages for each tool, as this will you give you the exact inputs that a tool requires.  
  
## Writing the workflow  
  
### Workflow template  
  We won't go into all the metadata you can provide to build a packaged bioinformatics tool, (you can find that out in the [building bioinformatics tools tutorial](#)), but here's a basic template:  
  
```python  
import janis as j  
import janis_bioinformatics  
  
# tool inputs go here  
  
class AlignSortedBam(j.Workflow):  
 def __init__(self):  
  super().__init__(self, "workflow_identifier_here")  
 # assemble our workflow here
```  

(Make sure to replace `"workflow_identifier_here"` with your own workflow name.)  
  
### Importing from `janis_bioinformatics`  
  
Throughout this workflow, we'll need to import `tools` and `data_types` from `janis_bioinformatics`. You can view the import path of a tool / data type in the documentation.  
  
We'd recommend you import from `janis_bioinformatics` instead of ~~`janis.bioinformatics`~~ due to the tricky way Python imports work.  
  
#### Importing Data Types  
  
We're going to need to import a series of `DataTypes` into our file, you can import `Bam`, `BamBai`, `Fastq`, `Sam`, `FastaWithDict` by placing the following line below the `"# tool inputs go here"` comment in the workflow template:   
  
```python  
from janis_bioinformatics.data_types import Bam, BamBai, Fastq, Sam, FastaWithDict  
```  
  
#### Importing tools  
  
> Pro tip: Import the exact version of the tool you want (over the latest) for maximum reproducibility!  
  
Let's import the tools that we've picked above by placing the following lines just under the data_type imports.  
```python  
from janis_bioinformatics.tools.bwa import BwaMem_0_7_15  
from janis_bioinformatics.tools.samtools import SamToolsView_1_7  
from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_0  
```  
  
#### Turning a `Tool` into a `Step`  
  
Now that we've imported our tools, within the `__init__` block of our subclassed workflow, we'll wrap all of our tools separately in a `janis.Step`.  
  
> (The `janis.Step` initialiser has the following contract):  
> ```python  
> janis.Step(identifier: str, tool: Tool, meta: Optional[Any]=None,    
 >  label: str = None, doc: str = None)  
> ```  
  
We'll give the step a name that corresponds with the tool we're going to use and assign the result to some variable. This makes editing our workflows easier, and also allows better tracking when we run the workflow.  
  
To keep things simple, we'll use the step identifier `"bwamem"` and the variable name `bwamem` for the Bwa mem step, which gives:
  
```python  
bwamem = j.Step("bwamem", BwaMem_0_7_15()) ```  
  
By repeating this process for the remainder of the tools, you should get the following code:  
  
```python  
def __init__(self):  
  super().__init__(self, "workflow_identifier_here")  
  
 # Steps bwa = j.Step("bwa", BwaMem_0_7_15())   samtools = j.Step("samtools", SamToolsView_1_7())    
   sortsam = j.Step("sortsam", Gatk4SortSam_4_0())  
```  
  
If so, nice work!  
  
### Declaring inputs  
  
By looking at the tool definitions, we can determine which inputs we need to provide the tools to get the workflow running.   
  
If you've used some of these tools before, you'll know that you may have had to provide an output name, fortunately `Janis` takes over that responsibility and generates a random filename.  
  
Let's analyse the inputs of our tools:

- BWA mem
    - Required
        - `reads` (`Fastq`)           
        - `reference` (`FastWithDict`)
    - Optional
        - `readGroupHeaderLine` (`String`)
- Samtools View
    - Required
        - `sam` (`Sam`)
- SortSam
    - Required:
        - `Bam` (`Bam`)         
        - `sortOrder` (`String`)
    - Optional
        - `createIndex` (`Boolean`)        
        - `validationStringency` (`String`)
        - `maxRecordsInRam` (`Int`)        
 

#### Putting this together  
  
We're going to connect some of the tools together, and default others, so we only need to expose the following three inputs:  
  
1. fastq (`Fastq`)  
2. reference (`FastaWithDict`)  
3. readGroupHeaderLine (`String`)  
  
The value in brackets is the DataType. For reference, the contract for `janis.Input` is the following.  
```python  
janis.Input(  
 identifier: str,            # unique identifier  
 data_type: DataType,         # data type (eg: String(), Int(), Bam()) value_from_input=None,       #    label: str=None,   
    doc: str=None)  
```  
  
The `identifier` should be descriptive and unique within the tool, and our `data_type` should reflect the input we're trying to mirror. Our Python variable name doesn't have to match the identifier, but it's recommended that they are or at least very similar. We can ignore the other options for now. These inputs can be placed below the step definitions within the `__init__` method:  
  
```python  
def __init__(...)  
 # step definitions read_group_header = j.Input("readGroupHeaderLine", String())   reference = j.Input("reference", FastaWithDict())    
   fastqs = j.Input("fastq", Fastq())  
```  
  
#### Outputs  
We'll specify the outputs after connecting the tools together.  
  
  
## Connecting the workflow  
  
Let's first look briefly how connections work in `janis`.  
  
### Overview  
  
> Learn more about edges: [Building connections](https://janis.readthedocs.io/en/latest/tutorials/buildingconnections.html)  
  
Now comes the connections! For this section we're going to use the following methods on the sub-workflow:  
  
- `add_edge(start, finish)`  
- `add_edges(*edges)`  
  
When constructing pipelines with python classes, we place an `Input` or `Step` in the `start` position, and either a `Step` or `Output` and `finish` position.   
  
We should full qualify our step nodes to ensure we're connecting the right attributes together, we do this with a `.` (period) followed by the attribute name, consider the following example:  
  
```python  
start: j.Input = Input("nameForStep", ...)  
step1: j.Step = Step(...)     # step1 has the 'name' && 'outp' attributes  
step2: j.Step = Step(...)     # step2 has the 'inp' attribute  
  
# Connect 'start' to 'step1' with the following line:  
self.add_edge(start, step1.name)  
self.add_edges([  
 (step1.outp, step2.inp)])  
```  
  
> Technical note: a `janis.Step` has a custom `__getattr__` method that validates the property is an input or output,  and then returns a format that 'add_edge' / 'add_edges' can interpret.  
  
To initially connect our tools, we're going to focus on connecting the inputs and steps (and worry about default values soon).  
  
### Connecting our bioinformatics tools  
  
#### Bwa mem  
  
We're going to connect the input `fastqs` to the input of `bwamem`, along with the reference file and a `readGroupHeaderLine`. The would result in the following connection:  
```python  
self.add_edges([    
    (fastqs, bwa.reads),    
    (read_group_header, bwa.readGroupHeaderLine),    
(reference, bwa.reference) ])  
```  
  
#### SamTools  
  
We want to connect the output of BWA Mem (`out`)  to the input of `bam`.  
  
```python  
self.add_edge(bwa.out, samtools.sam)  
```  
  
#### SortSam  
Similar to the previous, we connect to output of SamTools to the `bam` input of SortSam. Interestingly, the input of SortSam is `BamPair` (`.bam` + `.bai`), however we don't need to explicitly know that as `janis` will handle connecting secondary files where appropriate.  
  
```python  
self.add_edge(samtools.out, sortsam.bam)  
```  
  
#### Connected!  
  
Our workflow is now connected!  
  
  
### Cascading defaults  
  
Sometimes we want to provide default behaviour for the tool, in this case we'll specify an `Input` with a default value. This will still allow our user to provide a value, but can fall back if they choose they want the standard behaviour.   
  
```python  
w.add_edge(j.Input("myDefault", j.Int(), default=1), step.thatNeedsDefault)  
```  
  
#### Defaults for Workflow  
  
You can copy these defaults into your workflow:  
  
  
##### SortSam defaults  
  
```python  
self.add_edges([    
(j.Input("sortOrder", j.String(), default="coordinate"), sortsam.sortOrder),    
(j.Input("createIndex", j.Boolean(), default=True), sortsam.createIndex),    
(j.Input("validationStringency", j.String(), default="SILENT"), sortsam.validationStringency),    
(j.Input("maxRecordsInRam", j.Int(), default=5000000), sortsam.maxRecordsInRam), ])  
```  
  
### Outputs  
  
When connecting outputs, you only need to create an output node with a unique identifier. In this workflow, we only need to connect the output of SortSam to an output. Hence we have the following line at the end of our `__init__` method.  
  
```python  
self.add_edge(sortsam.out, Output("out"))  
```  
  
## Running your workflow  
  
> See the [running workflows](/tutorials/runningworkflows) guide for more information on exporting and running your workflows.  
  
### Adding input files  
  
To run your workflow, you'll need inputs! An easy way to do this as at the end of each required (non-default) input we've created, add a `, value="/path/to/file.ext"` after the data type. This will allow `janis` to create an input file for you which can be a little tedious when you have to specify every secondary file in WDL.  
  
```python  
reference = j.Input("reference", FastaWithDict(), value="/path/to/hg38/reference.fasta")  
```  
  This results in the following expert of WDL:  
 ```json{  
  "alignsortedbam.reference": "/path/to/hg38/reference.fasta",  
  "alignsortedbam.reference_amb": "/path/to/hg38/reference.fasta.amb",  
  "alignsortedbam.reference_ann": "/path/to/hg38/reference.fasta.ann",  
  "alignsortedbam.reference_bwt": "/path/to/hg38/reference.fasta.bwt",  
  "alignsortedbam.reference_dict": "/path/to/hg38/reference.dict",  
  "alignsortedbam.reference_fai": "/path/to/hg38/reference.fasta.fai",  
  "alignsortedbam.reference_pac": "/path/to/hg38/reference.fasta.pac",  
  "alignsortedbam.reference_sa": "/path/to/hg38/reference.fasta.sa",  
}  
```  
  
### Running your workflow with Cromwell  
You can run your workflow by first creating an instance and exporting to wdl (or cwl), we'll do this below the workflow in a `if __name__ == "__main__" block:  
```python  
if __name__ == "__main__":  
 w = AlignSortedBam() w.translate("wdl", to_disk=True, write_inputs_file=True)  
```  
  
This will automatically export the workflow to the `~/Desktop/alignsortedbam/wdl/` where you'll end up with:  
  
- `alignsortedbam.wdl`  
- `alignsortedbam-job.json`  
- `tools.zip`  
- `tools/`  
  
Where tools is a directory containing the `wdl` definitions of the tools we've used.  
  
Open a terminal to this directory, and run the following command to execute the workflow in Cromwell:  
> You may need to add `-Dconfig.file=/path/to/config.conf` if you're not using Dropbox  
```bash  
java -jar '/path/to/cromwell.jar' run 'alignsortedbam.wdl' -i 'alignsortedbam-job.json' -p 'tools.zip'  
```  
  
## Finished!  
  
Congratulations, you've just built and ran a a portable pipeline!  
  
  
#### Final result  
  
For reference, here's what your python file might look like:  
  
```python  
from janis import Workflow, Step, Output, Input, String, Boolean, Int  
from janis_bioinformatics.data_types import Fastq, FastaWithDict  
  
from janis_bioinformatics.tools.bwa import BwaMem_0_7_15  
from janis_bioinformatics.tools.samtools import SamToolsView_1_7  
from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_0  
  
class AlignSortedBam(Workflow):  
  
    def __init__(self):  
        super(AlignSortedBam, self).__init__("alignsortedbam")  
  
        bwa = Step("bwa", BwaMem_0_7_15()) 
        samtools = Step("samtools", SamToolsView_1_7()) 
        sortsam = Step("sortsam", Gatk4SortSam_4_0())  
        read_group_header = Input("readGroupHeaderLine", String()) 
        reference = Input("reference", FastaWithDict()) fastqs = Input("fastq", Fastq())         

        # S1: BWA mem  
        self.add_edges([  
            (fastqs, bwa.reads), 
            (read_group_header, bwa.readGroupHeaderLine), 
            (reference, bwa.reference) 
        ])  

        # S2: SamTools  
        self.add_edge(bwa.out, samtools.sam)  
  
        # S3: SortSam  
        self.add_edge(samtools.out, sortsam.bam)  

        # + S3 defaults  
        self.add_edges([  
            (Input("sortOrder", String(), default="coordinate"), sortsam.sortOrder),  
            (Input("createIndex", Boolean(), default=True), sortsam.createIndex),  
            (Input("validationStringency", String(), default="SILENT"), sortsam.validationStringency),  
            (Input("maxRecordsInRam", Int(), default=5000000), sortsam.maxRecordsInRam),  
        ])  
        # connect to output
        self.add_edge(sortsam.out, Output("out"))  

if __name__ == "__main__":  
    w = AlignSortedBam()
    w.translate("wdl", to_disk=True, write_inputs_file=True)  
 
```