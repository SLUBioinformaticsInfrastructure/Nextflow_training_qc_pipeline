#! /usr/bin/env nextflow

workflow{
	ch_input = Channel.fromFilePairs(params.input_read_pairs, checkIfExists: true )	

	// quality control
	FASTQC ( ch_input )

	// merging paired end reads
	FLASH2 ( ch_input )

	// read mapping to data base
	BOWTIE2 ( Channel.fromPath(params.bowtie2_db, checkIfExists: true ).toList(),
		 FLASH2.out.merged_reads) 

	// assembling reports
	MULTIQC ( FASTQC.out.qc_zip.collect(), FLASH2.out.log.collect(), BOWTIE2.out.logs.collect() )
}


// Publish directories are numbered to help understand processing order
// all variables named params.name are listed in params.yml

// fast quality control of fastq files
process FASTQC {

	input:
	tuple val(id), path(reads)


	// directives
	container 'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--hdfd78af_1' 
	publishDir "$params.outdir/01_fastqc" 

	script: 
	"""

	fastqc \\
		--noextract \\
		$reads

	"""

	output:
	path "${id}*fastqc.html", 	emit: qc_html
	path "${id}*fastqc.zip", 	emit: qc_zip

}


process FLASH2 {

	input: 
	tuple val(id), path(reads)

	// directives
	container 'https://depot.galaxyproject.org/singularity/flash2:2.2.00--hed695b0_2'
	publishDir "$params.outdir/02_flash"

	script: 
	"""
	
	flash2 \\
		$reads \\
		--output-prefix="${id}.flash2" \\
		--max-overlap=150 \\
		| tee -a ${id}_flash2.log


	"""

	output: 
	tuple val(id), path ("${id}.flash2.extendedFrags.fastq"), 	emit: merged_reads
	path "${id}.flash2.notCombined*.fastq", 	emit: notCombined_fastq
	path "${id}.flash2.hist", 			emit: flash2_hist
	path "${id}.flash2.histogram", 			emit: flash2_histogram
	path "${id}_flash2.log", 			emit: log
	
}


process BOWTIE2 {

	input: 
	path(bowtie2_db)
	tuple val(id), path(merged_reads)

	// directives:
	container 'https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h6fed5c7_2'
	publishDir "$params.outdir/03_bowtie2"

	script:
	db_name = bowtie2_db.find{it.name.endsWith(".1.bt2")}.name.minus(".1.bt2")
	"""

	bowtie2 \\
		-x $db_name \\
		-U $merged_reads \\
		-S ${id}_bowtie2_merged_${db_name}.sam \\
		--no-unal \\
		|& tee -a ${id}_bowtie2_merged_${db_name}.log

	"""

	output:
	path "${id}_bowtie2_merged_${db_name}.log", emit: logs
	path "${id}_bowtie2_merged_${db_name}.sam", emit: aligned_reads 
	

}

process MULTIQC {

	input: 
	path(fastqc_zips) 
	path(flash2_logs)
	path(bowtie2_log)

	// directives
	container 'https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0'
	publishDir "$params.outdir/04_multiqc"


	script: 
	"""
	multiqc \\
    		--force \\
    		--title "amrei test sub-set" \\
		.
	"""

	output: 
	path "*"
}
