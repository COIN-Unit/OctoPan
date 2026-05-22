version 1.0
import "./QualityControl.wdl" as QC
import "./MapToPangenome.wdl" as mapping
import "./SurjectBAM.wdl" as SurjectBAM
import "./VariantCalling.wdl" as Varcall
import "./Manta.wdl" as MantaVarcall
import "./STRipyPipeline.wdl" as STRipy
import "./SQUIRLS.wdl" as Squirls
import "./SNVstory.wdl" as snvstory
import "./ExtractChrM.wdl" as ExtractMitoChr
import "./SexDetermination.wdl" as SexDetermination
import "./MitoHPC.wdl" as MitoHPC
import "./Versions.wdl" as getVersions

workflow OctoPanWorkflow {
    input {
        File genderCheck
        File octopanToolkit_container
        File sample
        File json_to_csv_graph
        File samplesInfo
        File refGraph
        File haplo
        File SnpEff
        File grchDict
        File grchPath
        File chmPath
        File ref_hg38
        File ref_chm13
        File refAssembly
        File deepvariant_container
        File stripy_container
        File mitohpc_container
        File disease_gene_panel
        File squirls_Data
        File gnomAD
        File snv_resource
        File snvstory_container
        String mode
        Int cpu = 32
        Int mem_mb = 7100
        Int time = 1200
    }
    call QC.QualityControl as QualityCheck {
        input:
            octopanToolkit_container = octopanToolkit_container,
            sample = sample,
    }
    call mapping.MapToPangenomeWorkflow as PanMap {
        input:
            octopanToolkit_container = octopanToolkit_container,
            sample = sample,
            refGraph = refGraph,
            haplo = haplo,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }
    call SurjectBAM.SurjectBAMWorkflow as BamSurject {
        input:
            octopanToolkit_container = octopanToolkit_container,
            samplesInfo = samplesInfo,
            gafFile = PanMap.graph,
            refGraph = refGraph,
            chmPath = chmPath,
            grchPath = grchPath,
            grchDict = grchDict,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }
    call SexDetermination.SexDeterminationWorkflow as SexDet {
        input:
            octopanToolkit_container = octopanToolkit_container,
            genderCheck = genderCheck,
            bam_hg38 = bam_hg38,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }
    call Varcall.VariantCallingWorkflow as DeepVarCall {
        input:
            SnpEff = SnpEff,
            bam_hg38 = BamSurject.bam_hg38,
            bam_hg38_idx = BamSurject.bam_hg38_idx,
            bam_chm13 = BamSurject.bam_chm13,
            bam_chm13_idx = BamSurject.bam_chm13_idx,
            ref_hg38 = ref_hg38,
            ref_chm13 = ref_chm13,
            mode = mode,
            octopanToolkit_container = octopanToolkit_container,
            deepvariant_container = deepvariant_container,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }
    call MantaVarcall.MantaVariantCallingWorkflow as mantaVarCall {
            input:
                octopanToolkit_container = octopanToolkit_container,
                bam_hg38 = BamSurject.bam_hg38,
                bam_hg38_idx = BamSurject.bam_hg38_idx,
                bam_chm13 = BamSurject.bam_chm13,
                bam_chm13_idx = BamSurject.bam_chm13_idx,
                ref_hg38 = ref_hg38,
                ref_chm13 = ref_chm13,
                mode = mode,
                cpu = cpu,
                mem_mb = mem_mb,
                time = time
        }
    call STRipy.STRipyPipelineWorkflow as STRipyPipline {
        input:
            json_to_csv_graph = json_to_csv_graph,
            octopanToolkit_container = octopanToolkit_container,
            bam_hg38 = BamSurject.bam_hg38,
            bam_hg38_idx = BamSurject.bam_hg38_idx,
            refAssembly = refAssembly,
            samplesInfo = samplesInfo,
            stripy_container = stripy_container,
            sampleSex = SexDet.sampleSex,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }
    call Squirls.SQUIRLS as SquirlsPipline {
        input:
            SnpEff = SnpEff,
            octopanToolkit_container = octopanToolkit_container,
            vcf_file_hg38 = DeepVarCall.vcf_file_hg38,
            disease_gene_panel = disease_gene_panel,
            squirls_Data = squirls_Data,
            gnomAD = gnomAD,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }
    call snvstory.SNVstory as SNVstoryPipline {
        input:
            vcf_file_hg38 = DeepVarCall.vcf_file_hg38,
            snvstory_container = snvstory_container,
            snv_resource = snv_resource,
            mode = mode,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }
    call ExtractMitoChr.MitoSubset as subMito {
        input:
            bam_hg38 = BamSurject.bam_hg38,
            bam_hg38_idx = BamSurject.bam_hg38_idx,
            octopanToolkit_container = octopanToolkit_container,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    }

    call MitoHPC.MitoHPCWorkflow as MitoWorkflow {
        input:
            bam_hg38 = BamSurject.bam_hg38,
            bam_hg38_idx = BamSurject.bam_hg38_idx,
            mitohpc_container = mitohpc_container,
            octopanToolkit_container = octopanToolkit_container,
            cpu = cpu,
            mem_mb = mem_mb,
            time = time
    } 
    call getVersions.Versions as getV {
        input:
            SnpEff = SnpEff,
            squirls_Data = squirls_Data,
            octopanToolkit_container = octopanToolkit_container,
            deepvariant_container = deepvariant_container, 
            stripy_container = stripy_container,
            snvstory_container = snvstory_container,
            mem_mb = mem_mb,
            cpu = cpu,
            time = time
    }
    output { 
        File FastqcReportDir = QualityCheck.FastqcReportDir
        File MultiqcReportFile = QualityCheck.MultiqcReportFile
        File graph = PanMap.graph
        File bam_hg38 = BamSurject.bam_hg38
        File bam_hg38_idx = BamSurject.bam_hg38_idx
        File bam_chm13 = BamSurject.bam_chm13
        File bam_chm13_idx = BamSurject.bam_chm13_idx
        File sampleSex = SexDet.sampleSex
        File vcf_file_hg38 = DeepVarCall.vcf_file_hg38
        File vcf_file_chm13 = DeepVarCall.vcf_file_chm13
        File manta_vcf_file_hg38 = mantaVarCall.vcf_file_hg38
        File manta_vcf_file_chm13 = mantaVarCall.vcf_file_chm13
        File stripy_output = STRipyPipline.stripy_output
        File stripy_html = STRipyPipline.stripy_html
        File squirls_vcf = SquirlsPipline.squirls_vcf
        File squirls_html = SquirlsPipline.squirls_html        
        File snv_out = SNVstoryPipline.snv_out
        File chrm_bam = subMito.chrm_bam
        File toolsVersion = getV.toolsVersion
        File mtdnaVCF = MitoWorkflow.mtdnaVCF
        File haploReport = MitoWorkflow.haploReport
    }
} 
###############################################################################################################

 

 

 

 