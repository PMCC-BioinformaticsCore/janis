version 1.0

task RemoveAdapt {
  input{
    File sequence
    String adapter
    String sequence_name
  }

  Int disk_size = ceil(size(sequence, "GiB") * 2)
  Int memory_size = 6000

  command <<<
    cutadapt -a ~{adapter} -o ~{sequence_name}_trim.fastq.gz ~{sequence} --minimum-length 64
  >>>

  runtime {
    docker:"kfdrc/cutadapt"
    cpu:1
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "RemoveAdapt"
    mem:"~{memory_size}M"
    time:"10:00:00"
  }

  output {
    File trim_seq = "~{sequence_name}_trim.fastq.gz"
  }
}
