version 1.0

workflow MantaVariantCallingWorkflow {
    input {
        File octopanToolkit_container
        File bam_hg38
        File bam_hg38_idx
        File bam_chm13
        File bam_chm13_idx
        File ref_hg38
        File ref_chm13
        String mode
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000

    }
    call MantaVariantCallingHg38 {
        input:
            octopanToolkit_container = octopanToolkit_container,
            sample_id = sub(basename(bam_hg38), ".hprc-v1.1-mc-chm13_surject_hg38.bam$", ""),
            bam_hg38 = bam_hg38,
            bam_hg38_idx = bam_hg38_idx,
            ref_hg38 = ref_hg38,
            mem_mb = mem_mb,
            mode = mode,
            cpu = cpu,
            time = time
        }
    call MantaVariantCallingChm13 {
        input:
            octopanToolkit_container = octopanToolkit_container,
            sample_id = sub(basename(bam_chm13), ".hprc-v1.1-mc-chm13_surject_chm13.bam$", ""),
            bam_chm13 = bam_chm13,
            bam_chm13_idx = bam_chm13_idx,
            ref_chm13 = ref_chm13,
            mem_mb = mem_mb,
            mode = mode,
            cpu = cpu,
            time = time
        }
    output {
        File vcf_file_hg38 = MantaVariantCallingHg38.vcf_file_hg38
        File vcf_file_chm13 = MantaVariantCallingChm13.vcf_file_chm13
    }
}
##########################################################################################################
task MantaVariantCallingHg38 {
    input {
        String sample_id
        File bam_hg38
        File bam_hg38_idx
        File ref_hg38
        File octopanToolkit_container
        String mode
        Int cpu
        Int mem_mb
        Int time
    }
    command {
        Mode=~{mode}
        mkdir HG38_VCF
        singularity exec ~{octopanToolkit_container} samtools faidx ~{ref_hg38}
        # Step 1: Generate the configuration files
        singularity exec ~{octopanToolkit_container} configManta.py \
            --bam ~{bam_hg38} \
            --referenceFasta ~{ref_hg38} \
            --runDir HG38_VCF/~{sample_id}_manta \
            $( [[ "$Mode" == "WES" ]] && echo "--exome" )
        # Step 2: Run the Manta pipeline
        cd HG38_VCF/~{sample_id}_manta
        singularity exec ~{octopanToolkit_container} ./runWorkflow.py -m local

    }
    output {
        File vcf_file_hg38 = "HG38_VCF/~{sample_id}_manta/results/variants/diploidSV.vcf.gz"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
task MantaVariantCallingChm13 {
    input {
        String sample_id
        File bam_chm13
        File bam_chm13_idx
        File ref_chm13
        File octopanToolkit_container
        String mode
        Int cpu
        Int mem_mb
        Int time
    }
    command {
        Mode=~{mode}
        mkdir CHM13_VCF
        singularity exec ~{octopanToolkit_container} samtools faidx ~{ref_chm13}
        # Step 1: Generate the configuration files
        singularity exec ~{octopanToolkit_container} configManta.py \
            --bam ~{bam_chm13} \
            --referenceFasta ~{ref_chm13} \
            --runDir CHM13_VCF/~{sample_id}_manta \
            $( [[ "$Mode" == "WES" ]] && echo "--exome" )

        # Step 2: Run the Manta pipeline
        cd CHM13_VCF/~{sample_id}_manta
        singularity exec ~{octopanToolkit_container} ./runWorkflow.py -m local
    }
    output {
        File vcf_file_chm13 = "CHM13_VCF/~{sample_id}_manta/results/variants/diploidSV.vcf.gz"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}