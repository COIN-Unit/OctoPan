version 1.0

import "./IndexKmers.wdl" as Indexkmers


workflow MapToPangenomeWorkflow {
    input {
        File sample
        File refGraph
        File haplo
        File octopanToolkit_container
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000
    }

    call Indexkmers.IndexKmersWorkflow as KmerIndex {
        input:
            octopanToolkit_container = octopanToolkit_container,
            sample = sample,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }
    call MapToPangenome {
        input:
            sample = sample,
            sample_id = basename(sample),
            kff_file = KmerIndex.kffFile,
            RefGraph = refGraph,
            Haplotype = haplo,
            octopanToolkit_container = octopanToolkit_container,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }
    output {
            File graph = MapToPangenome.gaf_file
        }
}
###############################################################################################################
task MapToPangenome {
    input {
        File sample
        String sample_id
        File RefGraph
        File Haplotype
        File octopanToolkit_container
        File kff_file
        Int cpu
        Int mem_mb
        Int time
    }
    command <<<
    mkdir -p mapGraph
    fread=$(ls ~{sample}/*_1.fastq* 2>/dev/null | head -n 1)
    rread=$(ls ~{sample}/*_2.fastq* 2>/dev/null | head -n 1)
    singularity exec ~{octopanToolkit_container} \
      vg giraffe -Z ~{RefGraph} \
      -f "$fread" \
      -f "$rread" \
      --haplotype-name ~{Haplotype} \
      --kff-name ~{kff_file} \
      --index-basename ./ \
      -o gaf --sample ~{sample_id} --progress | singularity exec ~{octopanToolkit_container} bgzip > mapGraph/~{sample_id}.gaf.gz
    >>>
    output {
        File gaf_file = "mapGraph/~{sample_id}.gaf.gz"
    }

    runtime {
        cpu: cpu
        requested_memory_mb_per_core: mem_mb
        runtime_minutes: time
    }
}
