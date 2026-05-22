version 1.0

workflow Versions {
    input {
        File octopanToolkit_container
        File deepvariant_container
        File stripy_container
        File snvstory_container
        File SnpEff
        File squirls_Data
        Int cpu = 32
        Int mem_mb = 7000
        Int time = 4000
    }
    call versions {
        input:
            squirls_Data = squirls_Data,
            SnpEff = SnpEff,
            snvstory_container = snvstory_container,
            octopanToolkit_container = octopanToolkit_container,
            deepvariant_container = deepvariant_container,
            stripy_container = stripy_container,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    } 
    output {
        File toolsVersion = versions.toolsVersion
    }
}
###############################################################################################
task versions {
    input {
        File squirls_Data
        File SnpEff
        File snvstory_container
        File octopanToolkit_container
        File deepvariant_container
        File stripy_container
        Int cpu
        Int mem_mb
        Int time
    }
    command <<<
        echo "Fastqc -> " $(singularity exec ~{octopanToolkit_container} fastqc --version) > tools_version.txt
        echo "Multiqc -> " $(singularity exec ~{octopanToolkit_container} multiqc --version) >> tools_version.txt
        echo "Kmc -> " $(singularity exec ~{octopanToolkit_container} kmc | head -n 1) >> tools_version.txt
        echo "VG ->" $(singularity exec ~{octopanToolkit_container} vg 2>&1 | head -n 1) >> tools_version.txt
        echo "Samtools -> " $(singularity exec ~{octopanToolkit_container} samtools --version | head -n 1) >> tools_version.txt
        echo "Bgzip -> " $(singularity exec ~{octopanToolkit_container} bgzip -h | grep "Version") >> tools_version.txt
        echo "Tabix -> " $(singularity exec ~{octopanToolkit_container} tabix 2>&1 | grep "Version") >> tools_version.txt
        echo "DeepVariant ->" $(singularity exec ~{deepvariant_container} /opt/deepvariant/bin/run_deepvariant --version) >> tools_version.txt
        echo "SnpSift -> " $(java -jar ~{SnpEff}/SnpSift.jar 2>&1 | grep " version") >> tools_version.txt
        echo "Manta -> "  $(singularity exec ~{octopanToolkit_container} configManta.py | grep "Version") >> tools_version.txt
        echo "Stripy -> " $(singularity exec ~{stripy_container} cat /usr/local/bin/stripy-pipeline/stri.py | grep "version__ =" | cut -d '"' -f 2) >> tools_version.txt
        echo "Squirls -> " $(java -Xmx64G -jar ~{squirls_Data}/squirls-cli-2.0.1.jar -V) >> tools_version.txt
        echo "Snvstory -> " $(singularity run ~{snvstory_container} -h | grep "Ancestry Prediction") >> tools_version.txt
        echo "Haplocheck -> " $(singularity exec ~{octopanToolkit_container} haplocheck --version | head -n 2 | grep "h") >> tools_version.txt
        echo "MitoHPC -> 20240306"  >> tools_version.txt
    >>>

    output {
        File toolsVersion = "tools_version.txt"
    }

    runtime {
       cpu: cpu
       requested_memory_mb_per_core: mem_mb
       runtime_minutes: time
    }
}
