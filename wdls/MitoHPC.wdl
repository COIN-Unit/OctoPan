version 1.0

workflow MitoHPCWorkflow {
    input {
        File bam_hg38
        File bam_hg38_idx
        File mitohpc_container
        File octopanToolkit_container
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000

    } 
    call MitoHPC {
        input:
            bam_hg38 = bam_hg38,
            bam_hg38_idx = bam_hg38_idx,
            octopanToolkit_container = octopanToolkit_container,
            mitohpc_container = mitohpc_container,
            sample_id = sub(basename(bam_hg38), ".hprc-v1.1-mc-chm13_surject_hg38.bam$", ""),
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }

    call Haplocheck {
        input:
            octopanToolkit_container = octopanToolkit_container,
            mtdnaVCF = MitoHPC.mtdnaVCF,
            sample_id = sub(basename(bam_hg38), ".hprc-v1.1-mc-chm13_surject_hg38.bam$", ""),
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }
    output {
        File mtdnaVCF = MitoHPC.mtdnaVCF
        File haploReport = Haplocheck.haploReport
    }
}
##########################################################################################################
task MitoHPC {
    input {
        File octopanToolkit_container
        File bam_hg38
        File bam_hg38_idx
        File mitohpc_container
        String sample_id
        Int cpu
        Int mem_mb
        Int time
    }

    command <<<
        mkdir -p ~{sample_id}/out
        cp ~{bam_hg38} ./~{sample_id}/
        cp ~{bam_hg38_idx} ./~{sample_id}/

        singularity exec ~{octopanToolkit_container} samtools idxstats ./~{sample_id}/~{sample_id}".hprc-v1.1-mc-chm13_surject_hg38.bam" > ./~{sample_id}/~{sample_id}".hprc-v1.1-mc-chm13_surject_hg38.idxstats"
        singularity exec \
          --bind ./~{sample_id}:/mnt/bams \
          --bind ./~{sample_id}/out:/mnt/out \
          ~{mitohpc_container} \
          bash -c "
            export PATH=/MitoHPC/scripts:/MitoHPC/bin:/usr/local/bin:/usr/bin:/bin:$PATH
            export HP_SDIR=/MitoHPC/scripts
            export HP_ADIR=/mnt/bams
            export HP_ODIR=/mnt/out
            export HP_M=mutect2
            export HP_RNAME=hs38DH

            cd /MitoHPC/scripts
            source init.sh

            cd /mnt/bams
            bash /MitoHPC/scripts/filter.sh \
                ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38 \
                ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam \
                /mnt/out/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38
          "
    >>>

    output {
        File mtdnaVCF = "./~{sample_id}/out/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.mutect2.max.vcf"
    }

    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
task Haplocheck {
    input {
        File octopanToolkit_container
        File mtdnaVCF
        String sample_id
        Int cpu
        Int mem_mb
        Int time
    }
    command <<<
    singularity exec ~{octopanToolkit_container} haplocheck ~{mtdnaVCF} --out ./~{sample_id}".hprc-v1.1-mc-chm13_surject_hg38_haplo.haplocheck"

    >>>
    output {
        File haploReport = "./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_haplo.html"
    }
    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}  