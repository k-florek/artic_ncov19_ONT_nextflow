#!/usr/bin/env nextflow

//Description: Adaptation of ARTIC Network nCoV-2019 Bioinformatics SOP
//Available: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
//Author of this Nextflow: Kelsey Florek
//eMail: kelsey.florek@slh.wisc.edu

//starting parameters
params.fast5_dir = ""
params.run_prefix = "artic_ncov19"
params.maxcpus = 8


Channel
    .fromPath( "${params.fast5_dir}/**.fast5")
    .ifEmpty { exit 1, "Cannot find any fast5 files in: ${params.reads} Path must not end with /" }
    .set { raw_fast5 }

// Basecalling with Guppy
process guppy_basecalling {
  input:
    file(fast5s) from raw_fast5.collect()

  output:
    file "fastq/*.fastq" into fastq_reads
    file "fastq/sequencing_summary.txt" into seq_sum

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

// Demultiplexing and polishing with nanopolish
process artic_gather {
  input:
    file(reads) from fastq_reads.collect()
    file(summary) from seq_sum

  """
  artic gather --min-length 400 --max-length 700 --prefix ${params.run_name} --directory ./
  artic demultiplex --threads ${params.maxcpus} ${params.run_name}.fastq
  nanopolish index -s ${params.run_name}*sequencing_summary.txt -d ./ ${params.run_name}.fastq
  """
}
/*
process artic_demultiplex

}
*/
