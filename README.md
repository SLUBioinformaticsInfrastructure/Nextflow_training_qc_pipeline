# Nextflow_training_qc_pipeline
This is a short pipeline for students to improve their understanding of using Nextflow as a workflow manager. 



Tools used in the pipeline: 
- Fastqc for quality control
- Flash2 for merging paired end reads
- Bowtie2 for mapping the reads to a data base
- Multiqc to assemble the reports

Goals for training: 
- Show a simple pipeline that reproduces a workflow the students have already worked on.
- Talk about reproducibility and automation.
- Deepen their understanding of nextflow by letting them add their own process. 
- Talk about quality control and how to read the outputs of the used tools.


The file "nextflow_conda-env.yml" contains the recipe to create a shared environment for everyone. There are several ways to activate this environment: 
To create a uniform running environment make a shared conda environment:

*Use Conda:*

`conda env create --prefix /path/to/nextflow_conda-env.yml`

To activate this environment, use `conda activate /path/to/nextflow-env `.
To deactivate an active environment, use `conda deactivate`

 *Use Pixi:*

 `pixi init --import nextflow_conda-env.yml `

 To activate the environment, use `pixi shell`. 

 Note that during the pixi init command one could specify which paltform(s) the environment will be used in. The above command did not specify anything, so the resulting environment will be specific to the paltform it's created in. 
