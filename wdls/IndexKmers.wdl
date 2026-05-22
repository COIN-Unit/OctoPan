version 1.0
workflow IndexKmersWorkflow {
    input {
        File octopanToolkit_container
        File sample
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000
    }
    call Indexing{
     input:
       octopanToolkit_container = octopanToolkit_container,
       sample = sample,
       sample_id = basename(sample),
       mem_mb = mem_mb,
       cpu = cpu,
       time = time

    }
    output {
        File kffFile = Indexing.kmc_db
    }     
}
#############################################################################################
task Indexing {
  input{
    File octopanToolkit_container
    File sample
    String sample_id
    Int cpu
    Int mem_mb
    Int time
  }
  command <<<
  Spath=~{sample}
  sampleID=~{sample_id}
  for fread in "$Spath"/*_1.fastq*; do
    read=$(basename "$fread" | sed -E 's/_1\.fastq(\.gz)?//')  
    fread="$Spath/${read}_1.fastq"
    rread="$Spath/${read}_2.fastq"
    [[ -f "${fread}.gz" ]] && fread="${fread}.gz"
    [[ -f "${rread}.gz" ]] && rread="${rread}.gz"
    echo "$fread" >> "${sampleID}.reads.txt"
    echo "$rread" >> "${sampleID}.reads.txt"
  done
  singularity exec ~{octopanToolkit_container} kmc -k29 -m128 -t32 -okff @"${sampleID}.reads.txt" "${sampleID}" .
  >>>
  output {
    File kmc_db = "~{sample_id}.kff"
  }
  runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
   }
}