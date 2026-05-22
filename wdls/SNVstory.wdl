version 1.0

workflow SNVstory {
    input {
        File vcf_file_hg38
        File snv_resource
        File snvstory_container
        String mode
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000
    }
    call RunSNVstory {
        input:
            snvstory_container = snvstory_container,
            sample_id = sub(basename(vcf_file_hg38), ".hprc-v1.1-mc-chm13_surject_hg38_filter_PASS.dv.vcf.gz$", ""),
            vcf = vcf_file_hg38,
            snv_resource = snv_resource,
            mode = mode,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    } 
    output {
        File snv_out = RunSNVstory.snv_out
    }
}
###############################################################################################
task RunSNVstory {
    input {
        File snvstory_container
        String sample_id
        File snv_resource
        File vcf
        String mode 
        Int cpu
        Int mem_mb
        Int time

    }
    command <<<
            mkdir -p ~{sample_id}
            mkdir -p tmp_data

            singularity run \
                --bind ${PWD}/tmp_data:/data \
                ~{snvstory_container} \
                --path ~{vcf} \
                --resource ~{snv_resource} \
                --output-dir ~{sample_id} \
                --genome-ver 38 \
                --mode ~{mode}
        >>>

    output {
        File snv_out = "~{sample_id}/"
    }

    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
