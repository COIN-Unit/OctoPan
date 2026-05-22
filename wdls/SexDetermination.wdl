version 1.0

workflow SexDeterminationWorkflow {
    input {
        File genderCheck
        File octopanToolkit_container
        File bam_hg38
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000

    } 
    call SexDetermination {
        input:
            octopanToolkit_container = octopanToolkit_container,
            genderCheck = genderCheck,
            bam_hg38 = bam_hg38,
            sample_id = sub(basename(bam_hg38), ".bam$", ""),
            cpu = cpu,
            mem_mb = mem_mb,
            time = time,
    }
    output{
        File sampleSex = SexDetermination.sampleSex
    }

}
##########################################################################################################
task SexDetermination {
    input {
        File octopanToolkit_container
        File genderCheck
        File bam_hg38
        String sample_id
        Int cpu
        Int mem_mb
        Int time
    }
    command <<<
    singularity exec ~{octopanToolkit_container} python3 ~{genderCheck} ~{bam_hg38} > ~{sample_id}'.sex'
    >>>
    
    output {
        File sampleSex = "~{sample_id}.sex"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
