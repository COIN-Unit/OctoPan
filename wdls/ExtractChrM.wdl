version 1.0
workflow MitoSubset {
    input {
        File octopanToolkit_container
        File bam_hg38
        File bam_hg38_idx
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000
    }
    call ExtractChrM {
        input:
            octopanToolkit_container = octopanToolkit_container,
            sample_id = sub(basename(bam_hg38), ".hprc-v1.1-mc-chm13_surject_hg38.bam$", ""),
            bam = bam_hg38,
            bam_hg38_idx = bam_hg38_idx,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    } 
    output {
        File chrm_bam = ExtractChrM.chrm_bam
    }
}
###############################################################################################
task ExtractChrM {
    input {
        File octopanToolkit_container
        String sample_id
        File bam
        File bam_hg38_idx      
        Int cpu
        Int mem_mb
        Int time
    }
    command <<<
        singularity exec ~{octopanToolkit_container} samtools view -h -b ~{bam} chrM \
        > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_chrM.bam
    >>>

    output {
        File chrm_bam = "~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_chrM.bam"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
