version 1.0

workflow STRipyPipelineWorkflow {
    input {
        File octopanToolkit_container
        File bam_hg38
        File bam_hg38_idx
        File refAssembly
        File samplesInfo
        File stripy_container
        File json_to_csv_graph
        File sampleSex
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000
    }

    call STRipyPipeline {
        input:
            json_to_csv_graph = json_to_csv_graph,
            samplesInfo = samplesInfo,
            sampleSex = sampleSex,
            sample_id = sub(basename(bam_hg38), ".hprc-v1.1-mc-chm13_surject_hg38.bam$", ""),
            bam_hg38 = bam_hg38,
            bam_hg38_idx = bam_hg38_idx,
            octopanToolkit_container = octopanToolkit_container,
            refAssembly = refAssembly,
            stripy_container = stripy_container,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }
    
    output {
        File stripy_output = STRipyPipeline.stripy_output
        File stripy_html = STRipyPipeline.stripy_html
    }
}
##############################################################################################################
task STRipyPipeline {
    input {
        File json_to_csv_graph
        File octopanToolkit_container
        File samplesInfo
        File bam_hg38
        File bam_hg38_idx
        File refAssembly
        File sampleSex
        String sample_id
        File stripy_container
        Int cpu
        Int mem_mb
        Int time
    }

    command <<<
        #!/bin/bash
        singularity exec ~{octopanToolkit_container} samtools faidx ~{refAssembly}
        Sex=Sex=$(cat ~{sampleSex} | awk '{print $2}')

        singularity exec ~{stripy_container} python3 /usr/local/bin/stripy-pipeline/stri.py \
        --genome hg38 \
        --reference ~{refAssembly} \
        --output ./ --locus ABCD3,AFF2,AR,ARX_1,ARX_2,ATN1,ATXN1,ATXN10,ATXN2,ATXN3,ATXN7,ATXN8OS,BEAN1,\
        C9ORF72,CACNA1A,CBL,CNBP,COMP,DAB1,DIP2B,DMD,DMPK,FGF14,FMR1,FOXL2,FXN,GIPC1,GLS,HOXA13_1,\
        HOXA13_2,HOXA13_3,HOXD13,HTT,JPH3,LRP12,MARCHF6,NIPA1,NOP56,NOTCH2NLC,NUTM2B-AS1,PABPN1,PHOX2B,\
        PPP2R2B,PRDM12,RAPGEF2,RFC1,RILPL1,RUNX2,SAMD12,SOX3,STARD7,TBP,TBX1,TCF4,THAP11,TNRC6A,XYLT1,YEATS2,\
        ZIC2,ZIC3,CSTB,EIF4A3,PRNP,VWA1,ZFHX3 \
        --input ~{bam_hg38} \
        --sex ${Sex}

        cp /cbio/projects/003/saifeldeen/tools/NMDscan/containers/stripy-pipeline/json_to_csv_graph.R ./json_to_csv_graph.R

        sed -i 's/sample.stripy.json/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.json/g; s/sample.stripy.csv/~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.csv/g' \
        json_to_csv_graph.R

        Rscript json_to_csv_graph.R

        paste -d , <(cut -d, -f1 ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.csv) \
        <(cut -d, -f2- ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.csv | sed -E 's/([[:upper:]])/\L\1/g') \
        > ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.redcap.csv

        sed -i '1 s/\./_/g; 1 s/:/_/g' \
        ~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.redcap.csv

    >>>

    output {
        File stripy_output = "~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.redcap.csv"
        File stripy_html = "~{sample_id}.hprc-v1.1-mc-chm13_surject_hg38.bam.html"
    }

    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
