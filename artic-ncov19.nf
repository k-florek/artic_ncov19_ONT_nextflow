#!/usr/bin/env nextflow

//Description: Adaptation of ARTIC Network nCoV-2019 Bioinformatics SOP
//Available: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Authors of this Nextflow: Kelsey Florek and Abigail Shockey
//Email: kelsey.florek@slh.wisc.edu

// Input channels
Channel
    .value( "${params.primers}")
    .ifEmpty { exit 1, "Primers used must be included." }
    .into { polish_primers }

Channel
    .fromPath( "${params.fast5_dir}/*sequencing_summary*")
    .ifEmpty { exit 1, "Cannot find sequencing summary in: ${params.fast5_dir} Path must not end with /" }
    .set { sequencing_summary }

// If we have fast5 files then start with basecalling
if(params.basecalling){
  Channel
      .fromPath( "${params.fast5_dir}/**.fast5")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .into { raw_fast5; polish_fast5 }

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

// If we have already basecalled get fastq and fast5 for polishing
else {
  Channel
      .fromPath( "${params.fastq_dir}/*.fastq*")
      .ifEmpty { exit 1, "Cannot find any fastq files in: ${params.fastq_dir} Path must not end with /" }
      .set { fastq_reads }

  Channel
      .fromPath( "${params.fast5_dir}/**.fast5")
      .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.fast5_dir} Path must not end with /" }
      .set { polish_fast5 }
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

Channel
    .from(fastq_minion_pipeline)
    .combine(polish_fast5.collect())
    .combine(sequencing_summary)
    .set { polish_files }

if(params.polishing=="nanopolish"){
  process arctic_nanopolish_pipeline {
    publishDir "${params.outdir}/pipeline_nanopolish", mode: 'copy'

    input:
      val primers from nanopolish_primers
      tuple file(name), file(fast5s), file(sequencing_summary) from polish_files

    output:
      file "*{.primertrimmed.bam,.vcf,.variants.tab,.consensus.fasta}" into nanopolish_output

    script:
      """
      mkdir fast5s
      mv *.fast5 fast5s/

      filename=${name}
      samplename=\${filename%.*}

      artic minion --normalise ${params.normalise} --threads ${params.threadspipejob} --scheme-directory /artic-ncov2019/primer_schemes --fast5-directory ./fast5s --sequencing-summary ${sequencing_summary} --read-file ${name} nCoV-2019/${primers} \$samplename
      """
  }
}

else {
  process arctic_medaka_pipeline {
    publishDir "${params.outdir}/pipeline_nanopolish", mode: 'copy'

    input:
      val primers from nanopolish_primers
      tuple file(name), file(fast5s), file(sequencing_summary) from polish_files

    output:
      file "*{.primertrimmed.bam,.vcf,.variants.tab,.consensus.fasta}" into medaka_output

    script:
      """
      mkdir fast5s
      mv *.fast5 fast5s/

      filename=${name}
      samplename=\${filename%.*}

      artic minion --normalise ${params.normalise} --threads ${params.threadspipejob} --scheme-directory /artic-ncov2019/primer_schemes --fast5-directory ./fast5s --sequencing-summary ${sequencing_summary} --read-file ${name} nCoV-2019/${primers} \$samplename
      """
  }
}
