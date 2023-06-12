#!/bin/bash

#SBATCH -J test_nf_exon
#SBATCH --mem=4G
#SBATCH -o ./test_nf_exon_%J.out
#SBATCH -D /scratch/mblab/edwardkang/exon_nf/

# set up singularity
eval $(spack load --sh singularityce@3.8.0)
export SINGULARITY_CACHEDIR="/scratch/mblab/edwardkang/singularity/cache"

# set up nextflow environment
eval $(spack load --sh nextflow@22.10.4)

# run nextflow
nextflow run exonPipeline_array.nf -c conf/sampleData_array.config
