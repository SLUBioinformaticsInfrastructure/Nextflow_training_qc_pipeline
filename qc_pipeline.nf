#! /usr/bin/env nextflow

workflow{
//	ch_input = Channel.fromFilePairs(params.input_read_pairs, checkIfExists: true )	
	ch_input = channel.fromSRA(params.ids, apiKey:'01372e962e0d32291f4294ead0715461d608')

	// quality control
	FASTQC ( ch_input )

	// merging paired end reads
	FLASH2 ( ch_input )

	// read mapping to data base
	BOWTIE2 ( Channel.fromPath(params.bowtie2_db, checkIfExists: true ).toList(),
		 FLASH2.out.merged_reads ) 

	// 
	KRAKEN2 ( Channel.fromPath(params.kraken2_db, checkIfExists: true ).toList(),
		 ch_input ) 

	// 
	BRACKEN ( Channel.fromPath(params.kraken2_db, checkIfExists: true ).toList(),
		 KRAKEN2.out.kraken2_report ) 

	// 
	KRAKENTOOLS ( KRAKEN2.out.kraken2_report ) 

	// delete parts of the labels, to make the plots look nicer
	CLEANKRONA ( KRAKENTOOLS.out.converted,
		ch_input )

	// visualize the data in a krona plot
	KRONA ( CLEANKRONA.out.krona_cleaned ) 

	// Combine the krona plots into one single html file
	COMBINEKRONA ( CLEANKRONA.out.krona_cleaned.collect() )

	// assemble the reports
	MULTIQC ( FASTQC.out.qc_zip.collect(), FLASH2.out.log.collect(), 
		BOWTIE2.out.logs.collect(), KRAKEN2.out.logs.collect(),
		BRACKEN.out.logs.collect() )
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
	path "${id}.flash2.hist", 	emit: flash2_hist
	path "${id}.flash2.histogram", 	emit: flash2_histogram
	path "${id}_flash2.log", 	emit: log
}

process BOWTIE2 {

	input: 
	path(bowtie2_db)
	tuple val(id), path(merged_reads)

	// directives:
	container 'https://depot.galaxyproject.org/singularity/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
	publishDir "$params.outdir/03_bowtie2"

	script:
	db_name = bowtie2_db.find{it.name.endsWith(".rev.1.bt2")}.name.minus(".rev.1.bt2")
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

process KRAKEN2{

	input: 
	path(kraken2_db)
	tuple val(id), path(reads)

	// directives:
	container 'https://depot.galaxyproject.org/singularity/mulled-v2-8706a1dd73c6cc426e12dd4dd33a5e917b3989ae:c8cbdc8ff4101e6745f8ede6eb5261ef98bdaff4-0'
	publishDir "$params.outdir/04_kraken2"

	script:
	"""
	kraken2 \\
	    -db $kraken2_db \\
	    --paired \\
        --gzip-compressed \\
		--output ${id}_kraken2_taxonomy.txt \\
		--classified-out "${id}_kraken2_classified#.fq" \\
		--report ${id}_kraken2_report.txt \\
		$reads \\
	    |& tee -a ${id}_kraken2.log
	"""

    output:
	path "${id}_kraken2.log", emit: logs
	path "${id}_kraken2_classified*.fq", emit: kraken2_classified_reads 
	path "${id}_kraken2_taxonomy.txt", emit: kraken2_taxonomy
	path "${id}_kraken2_report.txt", emit: kraken2_report
}


process BRACKEN{

	input: 
	path(kraken2_db)
	path(kraken_report)

	// directives:
	container 'https://depot.galaxyproject.org/singularity/bracken:2.9--py38h2494328_0'
	publishDir "$params.outdir/05_bracken"

	script:
	id = kraken_report.name.minus("_kraken2_report.txt")
	"""
    bracken \\
        -d '${kraken2_db}' \\
        -i '${kraken_report}' \\
        -o '${id}_bracken_report.txt' \\
	    |& tee -a ${id}_bracken.log
	"""

    output:
    path "${id}_bracken_report.txt", emit: reports
	path "${id}_bracken.log", emit: logs	
}

process KRAKENTOOLS{

	input: 
	path(kraken_report)

	// directives:
	container 'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0'
	publishDir "$params.outdir/06_krakentools"

	script:
	id = kraken_report.name.minus("_kraken2_report.txt")
	"""
    kreport2krona.py \\
		-o '${id}_converted.txt' \\
		-r '${kraken_report}' \\
	    |& tee -a ${id}_kreport2krona.log
	"""

    output:
    path "${id}_converted.txt", emit: converted
	path "${id}_kreport2krona.log", emit: logs	
}

process CLEANKRONA {
	input: 
    path(kronareport)
	tuple val(id), path(reads)


	container 'https://depot.galaxyproject.org/singularity/krakentools:1.2--pyh5e36f6f_0'
	publishDir "$params.outdir/07_krakentools"

	script:
	"""
	sed -E 's/[a-z]__//g' $kronareport > ${id}_krona.txt
	"""

	output:
	path "${id}_krona.txt", emit: krona_cleaned
}



process KRONA{

	input: 
	path(converted)

	// directives:
	container 'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1'
	publishDir "$params.outdir/08_krona"

	script:
	id = converted.name.minus("_krona.txt")
	"""
    ktImportText \\
		-o '${id}_krona.html' \\
		'${converted}' \\
	    |& tee -a ${id}_krona.log
	"""

    output:
    path "${id}_krona.html", emit: reports
	path "${id}_krona.log", emit: logs	
}

process COMBINEKRONA {

	input: 
    path(kronareport)

	// directives:
	container 'https://depot.galaxyproject.org/singularity/krona:2.8.1--pl5321hdfd78af_1'
	publishDir "$params.outdir/09_compiledkrona"

	script:
	"""
    ktImportText \\
        $kronareport \\
        -o combined.krona.html

	"""

	output:
	path "combined.krona.html", emit: combined_krona_html
	
}

process MULTIQC {

	input: 
	path(fastqc_zips) 
	path(flash2_logs)
	path(bowtie2_log)
	path(kraken2_logs)
	path(bracken_logs)

	// directives
	container 'https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0'
	publishDir "$params.outdir/10_multiqc"


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
