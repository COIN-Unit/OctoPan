version 1.0

workflow SQUIRLS {
    input {
        File SnpEff
        File vcf_file_hg38
        File disease_gene_panel
        File octopanToolkit_container
        File gnomAD
        File squirls_Data
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000
    }

    call RunSQUIRLS {
        input:
            SnpEff = SnpEff,
            vcf = vcf_file_hg38,
            sample_id = sub(basename(vcf_file_hg38), ".hprc-v1.1-mc-chm13_surject_hg38_filter_PASS.dv.vcf.gz$", ""),
            disease_gene_panel = disease_gene_panel,
            octopanToolkit_container = octopanToolkit_container,
            squirls_Data = squirls_Data,
            gnomAD = gnomAD,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }
    output {
            File squirls_vcf = RunSQUIRLS.squirls_vcf
            File squirls_html = RunSQUIRLS.squirls_html        
        }
}
########################################################################################################
task RunSQUIRLS {
    input {
        File SnpEff
        File vcf
        File disease_gene_panel
        File octopanToolkit_container
        File gnomAD
        File squirls_Data
        String sample_id
        Int cpu
        Int mem_mb
        Int time

    }
    command <<<
        # step 7: run SQUIRLS
        singularity exec ~{octopanToolkit_container} tabix -p vcf ~{gnomAD}
        # filter for relevant genes (neuromuscular and neurodegenerative disorder superpanel)

        zcat ~{vcf} | java -Xmx64g -jar ~{SnpEff}/SnpSift.jar \
        intervals ~{disease_gene_panel} \
        > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes.dv.vcf
        # add SQUIRLS scores
        java -Xmx64G -jar ~{squirls_Data}/squirls-cli-2.0.1.jar annotate-vcf -d ~{squirls_Data} \
        ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes.dv.vcf ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls.dv --threads 32 -f vcf

        # filter variants (SQUIRLS score >= 0.5)

        cat ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls.dv.vcf | java -Xmx64g -jar ~{SnpEff}/SnpSift.jar filter "(SQUIRLS_SCORE >= 0.5)" \
        > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls_0.5.dv.vcf

        # add gnomAD v4.0 allele frequencies (exomes)

        java -Xmx64g -jar ~{SnpEff}/SnpSift.jar annotate -info AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_mid,AF_nfe,AF_sas,AF_remaining \
        ~{gnomAD} ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls_0.5.dv.vcf \
        > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls_0.5_gnomad_v4.dv.vcf

        # filter for rare variants (MAF<5%)

        cat ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls_0.5_gnomad_v4.dv.vcf | java -Xmx64g -jar ~{SnpEff}/SnpSift.jar filter \
        "(!(exists AF_afr) | (AF_afr < 0.05)) & \
        (!(exists AF_amr) | (AF_amr < 0.05)) & \
        (!(exists AF_asj) | (AF_asj < 0.05)) & \
        (!(exists AF_eas) | (AF_eas < 0.05)) & \
        (!(exists AF_fin) | (AF_fin < 0.05)) & \
        (!(exists AF_mid) | (AF_mid < 0.05)) & \
        (!(exists AF_nfe) | (AF_nfe < 0.05)) & \
        (!(exists AF_sas) | (AF_sas < 0.05)) & \
        (!(exists AF_remaining) | (AF_remaining < 0.05))" \
        > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls_0.5_gnomad_v4_rare_0.05.dv.vcf

        # create html squirls report for viewing

        java -Xmx64G -jar ~{squirls_Data}/squirls-cli-2.0.1.jar annotate-vcf -d ~{squirls_Data} \
        ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls_0.5_gnomad_v4_rare_0.05.dv.vcf ~{sample_id}_rare_0.05_squirls_score_0.5 -f html

        >>>
    output {
        File squirls_vcf = "~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_filter_PASS_nmd_genes_squirls_0.5_gnomad_v4_rare_0.05.dv.vcf"
        File squirls_html = "~{sample_id}_rare_0.05_squirls_score_0.5.html"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
