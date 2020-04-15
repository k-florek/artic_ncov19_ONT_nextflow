#!/usr/bin/env nextflow

//Description: Adaptation of ARTIC Network nCoV-2019 Bioinformatics SOP
//Available: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Authors of this Nextflow: Kelsey Florek and Abigail Shockey
//Email: kelsey.florek@slh.wisc.edu

//starting parameters
params.primers = ""
params.fast5_dir = ""
params.fastq_dir = ""
params.run_prefix = "artic_ncov19"
params.outdir = "artic_ncov19_results"
params.basecalling = "FALSE"

// If we have fast5 files then start with basecalling
if(params.basecalling){
  Channel
      .fromPath( "${params.fast5_dir}/**.fast5")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .into { raw_fast5; nanopolish_fast5; medaka_fast5 }

  Channel
      .fromPath( "${params.fast5_dir}/*sequencing_summary*")
      .ifEmpty { exit 1, "Cannot find sequencing summary in: ${params.fast5_dir} Path must not end with /" }
      .set { sequencing_summary }
  Channel
      .value( "${params.primers}")
      .ifEmpty { exit 1, "Primers used must be included." }
      .into { nanopolish_primers; medaka_primers}

  process guppy_basecalling {
    input:
      file(fast5s) from raw_fast5.collect()

    output:
      file "fastq/*.fastq" into fastq_reads

    script:
      if(params.basecalling_mode == "fast"){
        """
        guppy_basecaller -c /opt/ont/guppy/data/dna_r9.4.1_450bps_fast.cfg -i . -s fastq -x auto -r
        """
      }else{
        """
        guppy_basecaller -c /opt/ont/guppy/data/dna_r9.4.1_450bps_hac.cfg -i . -s fastq -x auto -r
        """
      }
  }
}

// If we have already basecalled get both fastq and fast5 files
else {
  Channel
      .fromPath( "${params.fastq_dir}/*.fastq*")
      .ifEmpty { exit 1, "Cannot find any fastq files in: ${params.fastq_dir} Path must not end with /" }
      .set { fastq_reads }

  Channel
      .fromPath( "${params.fast5_dir}/**.fast5")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .set { nanopolish_fast5 }

  Channel
      .fromPath( "${params.fast5_dir}/*sequencing_summary*")
      .ifEmpty { exit 1, "Cannot find sequencing summary in: ${params.fast5_dir} Path must not end with /" }
      .set { sequencing_summary }

  Channel
      .value( "${params.primers}")
      .ifEmpty { exit 1, "Primers used must be included." }
      .into { nanopolish_primers; medaka_primers}
}

process guppy_demultiplexing {
  publishDir "${params.outdir}/demultiplexing", mode: 'copy'

  input:
    file(fastqs) from fastq_reads.collect()

  output:
    path("output_directory/barcode*",type:'dir',maxDepth:1) into demultiplexed_reads

  script:
    """
      guppy_barcoder --require_barcodes_both_ends -i . -s output_directory --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
    """
}

process artic_guppyplex {
  publishDir "${params.outdir}/guppyplex", mode: 'copy'
  stageOutMode = 'link'

  input:
    file(reads) from demultiplexed_reads.flatten()

  output:
    file "${params.run_prefix}_barcode*.fastq" into fastq_minion_pipeline

  script:
    """
    artic guppyplex \
    --min-length ${params.min_length} \
    --max-length ${params.max_length} \
    --directory barcode* \
    --prefix ${params.run_prefix}
    """
}

process arctic_minion_pipeline {
  publishDir "${params.outdir}/pipeline_nanopolish", mode: 'copy'
  stageInMode = 'link'

  input:
    val primers from nanopolish_primers
    file(name) from fastq_minion_pipeline
    file(fast5s) from nanopolish_fast5.collect()
    file(sequencing_summary) from sequencing_summary

  output:
    file "*{.primertrimmed.bam,.vcf,.variants.tab,.consensus.fasta}" into output 
  script:
    """
    mkdir fast5s
    mv *.fast5 fast5s/

    filename=${name}
    samplename=\${filename%.*}

    artic minion --normalise ${params.normalise} --threads ${params.threadspipejob} --scheme-directory /artic-ncov2019/primer_schemes --fast5-directory ./fast5s --sequencing-summary ${sequencing_summary} --read-file ${name} nCoV-2019/${primers} \$samplename
    """
}
