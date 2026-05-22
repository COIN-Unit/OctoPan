version 1.0
workflow QualityControl {
  input{
    File octopanToolkit_container
    File sample
 }
 call Fastqc{
  input: 
    octopanToolkit_container = octopanToolkit_container,
    sample = sample
  }
 call Multiqc{
  input:
    octopanToolkit_container = octopanToolkit_container,
    FastqcOUt = Fastqc.fastqcReports
  }
  output {
    File FastqcReportDir = Fastqc.fastqcReports
    File MultiqcReportFile = Multiqc.multiqcOut
  }
}
##############################################################
task Fastqc {
  input{
    File sample
    File octopanToolkit_container
  }
  command <<<
  Spath=~{sample}
  sampleID=$(basename "$Spath")
  mkdir -p FASTQC_Reports/"${sampleID}_QC"
  for fread in "$Spath"/*_1.fastq*; do
    read=$(basename "$fread" | sed -E 's/_1\.fastq(\.gz)?//')
    if [[ "$fread" == *.gz ]]; then
      rread="$Spath/${read}_2.fastq.gz"
    else
      rread="$Spath/${read}_2.fastq"
    fi
    fread="$Spath/${read}_1.fastq"
    [[ -f "$fread.gz" ]] && fread="$fread.gz"
    [[ -f "$rread.gz" ]] && rread="$rread.gz"
    if [[ -f "$fread" && -f "$rread" ]]; then
      singularity exec ~{octopanToolkit_container} fastqc "$fread" "$rread" --outdir FASTQC_Reports/"${sampleID}_QC" --threads 32
    else
      echo "Skipping $read: missing one of the pair files"
    fi
  done
  >>>
  output {
    File fastqcReports = "FASTQC_Reports"
  }
}
task Multiqc {
    input{
        File octopanToolkit_container
        File FastqcOUt
    }
    command <<<
    mkdir "Multiqc_Report"
    singularity exec ~{octopanToolkit_container} multiqc -o Multiqc_Report ~{FastqcOUt}"/."
    >>>
    output {
        File multiqcOut = "Multiqc_Report"
    }
}