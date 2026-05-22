version 1.0

workflow VariantCallingWorkflow {
    input {
        File octopanToolkit_container
        File bam_hg38
        File bam_hg38_idx
        File bam_chm13
        File bam_chm13_idx
        File ref_hg38
        File ref_chm13
        File SnpEff
        File deepvariant_container
        String mode
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000

    }
    call VariantCallingHg38 {
        input:
            sample_id = sub(basename(bam_hg38), ".hprc-v1.1-mc-chm13_surject_hg38.bam$", ""),
            bam_hg38 = bam_hg38,
            bam_hg38_idx = bam_hg38_idx,
            ref_hg38 = ref_hg38,
            mode = mode,
            SnpEff = SnpEff,
            octopanToolkit_container = octopanToolkit_container,
            deepvariant_container = deepvariant_container,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }

    call VariantCallingChm13 {
        input:
            octopanToolkit_container = octopanToolkit_container,
            sample_id = sub(basename(bam_chm13), ".hprc-v1.1-mc-chm13_surject_chm13.bam$", ""),
            bam_chm13 = bam_chm13,
            bam_chm13_idx = bam_chm13_idx,
            ref_chm13 = ref_chm13,
            SnpEff = SnpEff,
            mode = mode,
            deepvariant_container = deepvariant_container,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }
    output {
        File vcf_file_hg38 = VariantCallingHg38.vcf_file_hg38
        File vcf_file_chm13 = VariantCallingChm13.vcf_file_chm13
    }
}
##########################################################################################################
task VariantCallingHg38 {
    input {
        File octopanToolkit_container
        String sample_id
        File bam_hg38
        File SnpEff
        File bam_hg38_idx
        File ref_hg38
        File deepvariant_container
        String mode
        Int cpu
        Int mem_mb
        Int time
    }
    command {
        mkdir HG38_VCF
        singularity exec ~{octopanToolkit_container} samtools faidx ~{ref_hg38}
        singularity exec ~{deepvariant_container} \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=~{mode} \
        --ref=~{ref_hg38} \
        --reads=~{bam_hg38} \
        --output_vcf=HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.dv.vcf.gz \
        --output_gvcf=HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.dv.g.vcf.gz \
        --sample_name=~{sample_id} \
        --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true" --num_shards=32

    
        zcat HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.dv.vcf.gz | java -Xmx64g -jar ~{SnpEff}/SnpSift.jar filter "(FILTER = 'PASS')" \
        > HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS.dv.vcf

        singularity exec ~{octopanToolkit_container} bgzip -c HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS.dv.vcf \
        > HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS.dv.vcf.gz
        
        singularity exec ~{octopanToolkit_container} tabix -p vcf HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS.dv.vcf.gz
    }
    output {
        File vcf_file_hg38 = "HG38_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS.dv.vcf.gz"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
task VariantCallingChm13 {
    input {
        String sample_id
        File octopanToolkit_container
        File bam_chm13
        File bam_chm13_idx
        File ref_chm13
        File SnpEff
        File deepvariant_container
        String mode
        Int cpu
        Int mem_mb
        Int time
    }
    command {
        mkdir CHM13_VCF
        singularity exec ~{octopanToolkit_container} samtools faidx ~{ref_chm13}
        singularity exec ~{deepvariant_container} \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=~{mode} \
        --ref=~{ref_chm13} \
        --reads=~{bam_chm13} \
        --output_vcf=CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.dv.vcf.gz \
        --output_gvcf=CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.dv.g.vcf.gz \
        --sample_name=~{sample_id} \
        --make_examples_extra_args="min_mapping_quality=1,keep_legacy_allele_counter_behavior=true,normalize_reads=true" --num_shards=32


        zcat CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.dv.vcf.gz | java -Xmx64g -jar ~{SnpEff}/SnpSift.jar filter "(FILTER = 'PASS')" \
        > CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_filter_PASS.dv.vcf

        singularity exec ~{octopanToolkit_container} bgzip -c CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_filter_PASS.dv.vcf \
        > CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_filter_PASS.dv.vcf.gz
        
        singularity exec ~{octopanToolkit_container} tabix -p vcf CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_filter_PASS.dv.vcf.gz
    }
    output {
        File vcf_file_chm13 = "CHM13_VCF/~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_filter_PASS.dv.vcf.gz"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
