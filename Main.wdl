version 1.0
import "./wdls/QualityControl.wdl" as QC
import "./wdls/MapToPangenome.wdl" as mapping
import "./wdls/SurjectBAM.wdl" as SurjectBAM
import "./wdls/SexDetermination.wdl" as SexDetermination
import "./wdls/VariantCalling.wdl" as Varcall
import "./wdls/Manta.wdl" as MantaVarcall
import "./wdls/STRipyPipeline.wdl" as STRipy
import "./wdls/SQUIRLS.wdl" as Squirls
import "./wdls/SNVstory.wdl" as snvstory
import "./wdls/ExtractChrM.wdl" as ExtractMitoChr
import "./wdls/MitoHPC.wdl" as MitoHPC
import "./wdls/Versions.wdl" as getVersions
import "./wdls/OctoPan.wdl" as OctoPan

workflow Main {
    input {
        String wd
        File samplesList
        File octopanToolkit_container
        File samplesInfo
        File refGraph
        File haplo
        File grchDict
        File grchPath
        File chmPath
        File ref_hg38
        File ref_chm13
        File refAssembly
        File SnpEff
        File deepvariant_container
        File stripy_container
        File disease_gene_panel
        File snvstory_container
        File mitohpc_container
        File squirls_Data
        File gnomAD
        File snv_resource
        File genderCheck
        File json_to_csv_graph
        String mode
        Int cpu = 32
        Int mem_mb = 7100
        Int time = 1200
    }
    call getSample{
        input:
            workingDir = wd,
            sampleFile = samplesList
    }
    scatter (sample in getSample.samplesPathFiles){
        call OctoPan.OctoPanWorkflow as MainOctoPan {
            input:
                sample = sample,
                samplesInfo = samplesInfo,
                refGraph = refGraph,
                haplo = haplo,
                grchDict = grchDict,
                grchPath = grchPath,
                chmPath = chmPath,
                ref_hg38 = ref_hg38,
                ref_chm13 = ref_chm13,
                json_to_csv_graph = json_to_csv_graph,
                refAssembly = refAssembly,
                mitohpc_container = mitohpc_container,
                snvstory_container = snvstory_container,
                octopanToolkit_container = octopanToolkit_container,
                deepvariant_container = deepvariant_container,
                stripy_container = stripy_container,
                disease_gene_panel = disease_gene_panel,
                squirls_Data = squirls_Data,
                gnomAD = gnomAD,
                snv_resource = snv_resource,
                genderCheck = genderCheck,
                SnpEff = SnpEff,
                mode = mode,
                cpu = cpu,
                mem_mb = mem_mb,
                time = time
        }
    }
output {
    Array[File] FastqcReportDir = MainOctoPan.FastqcReportDir
    Array[File] MultiqcReportFile = MainOctoPan.MultiqcReportFile
    Array[File] graph = MainOctoPan.graph
    Array[File] bam_hg38 = MainOctoPan.bam_hg38
    Array[File] bam_chm13 = MainOctoPan.bam_chm13
    Array[File] bam_hg38_idx = MainOctoPan.bam_hg38_idx
    Array[File] bam_chm13_idx = MainOctoPan.bam_chm13_idx
    Array[File] vcf_file_hg38 = MainOctoPan.vcf_file_hg38
    Array[File] vcf_file_chm13 = MainOctoPan.vcf_file_chm13
    Array[File] manta_vcf_file_hg38 = MainOctoPan.manta_vcf_file_hg38
    Array[File] manta_vcf_file_chm13 = MainOctoPan.manta_vcf_file_chm13
    Array[File] stripy_output = MainOctoPan.stripy_output
    Array[File] stripy_html = MainOctoPan.stripy_html
    Array[File] squirls_vcf = MainOctoPan.squirls_vcf
    Array[File] squirls_html = MainOctoPan.squirls_html        
    Array[File] snv_out = MainOctoPan.snv_out
    Array[File] chrm_bam = MainOctoPan.chrm_bam
    Array[File] mtdnaVCF = MainOctoPan.mtdnaVCF
    Array[File] haploReport = MainOctoPan.haploReport
    Array[File] toolsVersion = MainOctoPan.toolsVersion
    Array[File] sampleSex = MainOctoPan.sampleSex
  }
}
#########################################################################################
task getSample{
  input{
    String workingDir
    File sampleFile
  }
  command <<<
  for sample in $(cat ~{sampleFile}); do
    echo ~{workingDir}"/raw_reads/"${sample} | sort
  done > samplesPath.txt
  >>>
  output{
    Array[File] samplesPathFiles = read_lines("samplesPath.txt")
  }
}