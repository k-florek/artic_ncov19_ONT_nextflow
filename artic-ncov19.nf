#!/usr/bin/env nextflow

//Description: Adaptation of ARTIC Network nCoV-2019 Bioinformatics SOP
//Available: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Author of this Nextflow: Kelsey Florek
//eMail: kelsey.florek@slh.wisc.edu

//starting parameters
params.fast5_dir = ""
params.run_prefix = "artic_ncov19"
params.outdir = "artic_ncov19_results"


Channel
    .fromPath( "${params.fast5_dir}/**.fast5")
    .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.reads} Path must not end with /" }
    .into { raw_fast5; polish_fast5 }

process guppy_basecalling {
  input:
    file(fast5s) from raw_fast5.collect()

  output:
    file "fastq/*.fastq" into fastq_reads

  script:
  if(params.basecalling_mode == "fast"){
    """
    guppy_basecaller -c /opt/ont/guppy/data/dna_r9.4.1_450bps_fast.cfg -i ./ -s fastq -x auto -r
    """
  }else{
    """
    guppy_basecaller -c /opt/ont/guppy/data/dna_r9.4.1_450bps_hac.cfg -i ./ -s fastq -x auto -r
    """
  }
}

process artic_gather {
  input:
    file(reads) from fastq_reads.collect()

  output:
    file("${params.run_prefix}_fastq_pass.fastq") into fastq_polishing, fastq_demultiplexing
    file("${params.run_prefix}_sequencing_summary.txt") into summary

  script:
  """
  artic gather \
  --min-length ${params.min_length} \
  --max-length  ${params.max_length} \
  --prefix ${params.run_prefix} \
  --directory ./
  """
}

process artic_demultiplex {
  input:
    file(fastq_reads) from fastq_demultiplexing

  output:
    file("${params.run_prefix}_pass_*.fastq") into demultiplexed_reads

  script:
  """
  artic demultiplex \
  --threads ${params.threadsmultiplexjob} \
  ${fastq_demultiplexing}
  """
}

process artic_nanopolish {
  input:
    file(fastq_passing) from fastq_polishing
    file(seq_summary) from summary
    path fast5s, stageAs:'fast5/*' from polish_fast5.collect()

  output:
    file "*.index*" into nanopolish_indexs
    file "${params.run_prefix}_pass.fastq" into nano_passed_reads

  script:
  """
  nanopolish index \
  -s ${seq_summary} \
  -d fast5/
  ${fastq_passing}
  """
}

process artic_pipeline {
  publishDir "${params.outdir}", mode: "copy"

  input:
    file(nano_index) from nanopolish_indexs.collect()
    file(read_file) from demultiplexed_reads
    file(nano_reads) from nano_passed_reads

  output:
    file "*{.primertrimmed.bam,.vcf,.variants.tab,.consensus.fasta}" into output

  script:
  """
  filename=${read_file} && tmp=\${filename#*pass_} && samplename=\${tmp%.*}
  artic minion \
  --normalise ${params.normalise} \
  --threads ${params.threadspipejob} \
  --scheme-directory /artic-ncov2019/primer_schemes \
  --read-file ${read_file} \
  --nanopolish-read-file ${nano_reads} \
  nCoV-2019 \
  \$samplename
  """
}
