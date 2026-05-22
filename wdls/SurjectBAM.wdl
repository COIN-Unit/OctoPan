version 1.0

workflow SurjectBAMWorkflow {
    input {
        File samplesInfo
        File gafFile
        File refGraph
        File chmPath
        File grchPath
        File grchDict
        File octopanToolkit_container                
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000

    }
    call SurjectBAM {
        input:
            samplesInfo = samplesInfo,
            refGraph = refGraph,
            chmPath = chmPath,
            grchPath = grchPath,
            grchDict = grchDict,
            sample_id = sub(basename(gafFile), ".gaf.gz$", ""),
            gaf_file = gafFile,
            octopanToolkit_container = octopanToolkit_container,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
        }

    output {
        File bam_hg38 = SurjectBAM.bam_hg38
        File bam_hg38_idx = SurjectBAM.bam_hg38_idx
        File bam_chm13 = SurjectBAM.bam_chm13
        File bam_chm13_idx = SurjectBAM.bam_chm13_idx
    }


}
#############################################################################################################
task SurjectBAM {
    input {
        File samplesInfo
        File gaf_file
        File refGraph
        File chmPath
        File grchPath
        File grchDict
        String sample_id
        File octopanToolkit_container
        Int cpu
        Int mem_mb
        Int time
    }

    command <<<
        flowCellID=$(grep "^~{sample_id}" ~{samplesInfo} | awk '{print $3}')
        Lane=$(grep "^~{sample_id}" ~{samplesInfo} | awk '{print $4}')
        Lib=$(grep "^~{sample_id}" ~{samplesInfo} | awk '{print $5}')
        PLATFORM=$(grep "^~{sample_id}" ~{samplesInfo} | awk '{print $6}')
        singularity exec ~{octopanToolkit_container} \
        vg surject -x ~{refGraph} -G ~{gaf_file} --interleaved -F ~{grchPath} -b -N ~{sample_id} \
        -R 'ID:$flowCellID.$Lane LB:$Lib SM:~{sample_id} PL:$PLATFORM' > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam
        
        singularity exec ~{octopanToolkit_container} \
        vg surject -x ~{refGraph} -G ~{gaf_file} --interleaved -F ~{chmPath} -b -N ~{sample_id} \
        -R 'ID:$flowCellID.$Lane LB:$Lib SM:~{sample_id} PL:$PLATFORM' > ~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.bam

        #fix the hg38 BAM header
        singularity exec ~{octopanToolkit_container} \
        samtools reheader ~{grchDict} ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam > ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_fix.bam
        mv ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_fix.bam  ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam

        singularity exec ~{octopanToolkit_container} \
        samtools sort ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam -O BAM -o ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_sort.bam --threads 8
        
        singularity exec ~{octopanToolkit_container} \
        samtools index ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_sort.bam -@ 8

        singularity exec ~{octopanToolkit_container} \
        samtools sort ./~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.bam -O BAM -o ./~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_sort.bam --threads 8

        singularity exec ~{octopanToolkit_container} \
        samtools index ./~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_sort.bam -@ 8

        # change the chromosome names so that they are compatable with other tools
        singularity exec ~{octopanToolkit_container} \
        samtools view -H ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_sort.bam > ~{sample_id}.surject_hg38_bam_header.txt
        
        sed -i 's/GRCh38#0#//g' ~{sample_id}.surject_hg38_bam_header.txt
        
        singularity exec ~{octopanToolkit_container} \
        samtools reheader ~{sample_id}.surject_hg38_bam_header.txt ./~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_sort.bam > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam

        singularity exec ~{octopanToolkit_container} \
        samtools index ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam

        singularity exec ~{octopanToolkit_container} \
        samtools view -H ./~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_sort.bam > ~{sample_id}.surject_chm13_bam_header.txt

        singularity exec ~{octopanToolkit_container} \
        sed -i 's/CHM13#0#//g' ~{sample_id}.surject_chm13_bam_header.txt

        singularity exec ~{octopanToolkit_container} \
        samtools reheader ~{sample_id}.surject_chm13_bam_header.txt ./~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_sort.bam > ~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.bam

        singularity exec ~{octopanToolkit_container} \
        samtools index ~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.bam

        rm ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38_sort*
        rm ~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13_sort*
        >>>

    output {
        File bam_hg38 = "~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam"
        File bam_hg38_idx = "~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.bai"
        File bam_chm13 = "~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.bam"
        File bam_chm13_idx = "~{sample_id}.hprc-v1.1-mc-chm13_surject_chm13.bam.bai"
    }

    runtime {
           cpu: cpu
           requested_memory_mb_per_core: mem_mb
           runtime_minutes: time
    }
}

