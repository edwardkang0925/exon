params.trait = 'BMI'
params.numSlice = 2
params.minBinWidth = 45
params.minCPM = 3
params.minPercentsExpressed = 0.015
params.number_of_plates = 28
params.visitcode = "v1"
params.binLevelPval = 0.01

nextflow.enable.dsl=2

process PreprocessData {
    container './exon_latest.sif'

    output:
    path("outputs/preprocessed/pc_${params.trait}.csv")
    path("outputs/preprocessed/vst_${params.trait}.RDS")

    script:
    """
    Rscript /project/scripts/preprocess.R --trait ${params.trait} --minBinWidth ${params.minBinWidth} \
        --minCPM ${params.minCPM} --minPercentsExpressed ${params.minPercentsExpressed} \
        --binFile "/project/data/bins_sample.RDS"
    """
}

process PrepareRegression {
    container './exon_latest.sif'

    input:
    path(pc)
    path(vst)

    output:
    path("outputs/beforeRegression/${params.trait}/sliced/exonSlice_*.RDS")
    path("outputs/beforeRegression/${params.trait}/phenotype_df_beforeRegression_${params.trait}.RDS")
    path("outputs/beforeRegression/${params.trait}/combined_df_beforeRegression_${params.trait}.RDS")

    script:
    """
    Rscript /project/scripts/prepareRegression.R --trait ${params.trait} --numSlice ${params.numSlice} \
        --number_of_plates ${params.number_of_plates} --visitcode ${params.visitcode} \
        --pcFilePath ${pc} --vstFilePath ${vst}
    """
}

process StepwiseRegression {
    container './exon_latest.sif'
    label 'SLURM_array'

    input:
    path(slicedBinCount)
    path(phenotypeBeforeRegression)
    path(combinedDF)

    output:
    path("outputs/afterRegression/${params.trait}/sliced/")

    script:
    """
    Rscript /project/scripts/stepwise_regression_perbin.R \
        --combined_df ${combinedDF} \
        --exon_expression_path ${slicedBinCount} \
        --output_directory "outputs/afterRegression/${params.trait}/sliced/"
    """
}

process MergeResiduals {
    container './exon_latest.sif'
    publishDir "./", mode: 'copy'

    input:
    path(residualBinExpressionDir)

    output:
    path("outputs/afterRegression/${params.trait}/combined/residualExonExpression.RDS")

    script:
    """
    Rscript /project/scripts/combine_residuals.R \
        --input_prefix ${residualBinExpressionDir} \
        --output_prefix "outputs/afterRegression/${params.trait}/combined/"
    """
}

process PrepareConditionalGenesis {
    container './exon_latest.sif'

    input:
    path(mergedResidualsExpression)

    output:
    path("outputs/conditionalGenesisInput/${params.trait}/conditionalSlice_*.RDS")

    script:
    """
    Rscript /project/scripts/preprocessForConditional_all.R \
        --trait ${params.trait} \
        --numSlice ${params.numSlice} \
        --pvalThreshold ${params.binLevelPval} \
        --mergedResidualsPath ${mergedResidualsExpression} \
        --outputDir "outputs/conditionalGenesisInput/${params.trait}/"
    """
}

process RunConditionalGenesis {
    container './exon_latest.sif'
    label 'SLURM_array'

    input:
    path(conditionalGenesisInputSliceFile)

    output:
    path("outputs/conditionalGenesisOutput/${params.trait}/")

    script:
    """
    Rscript /project/scripts/conditionalGenesis_chunk.R \
        --trait ${params.trait} \
        --exon.path ${conditionalGenesisInputSliceFile} \
        --outputDir "outputs/conditionalGenesisOutput/${params.trait}/"
    """
}

process MergeConditionalGenesis {
    container './exon_latest.sif'
    publishDir "./", mode: 'copy'

    input:
    path(conditionalGenesisOutputDir)

    output:
    path("outputs/combinedConditionalGenesis/${params.trait}/comnbinedConditionalGenesis_${params.trait}.RDS")

    script:
    """
    Rscript /project/scripts/combineConditionalGenesis.R \
        --trait ${params.trait} \
        --inputDir ${conditionalGenesisOutputDir} \
        --outputFile "outputs/combinedConditionalGenesis/${params.trait}/comnbinedConditionalGenesis_${params.trait}.RDS"
    """
}


workflow {
    preprocessDataOut = PreprocessData()
    prepareRegressionOut = PrepareRegression(preprocessDataOut[0], preprocessDataOut[1])
    stepwiseRegressionOut = StepwiseRegression(prepareRegressionOut[0], prepareRegressionOut[1], prepareRegressionOut[2])

    // Collect all outputs from StepwiseRegression
    stepwiseRegressionOut.collect().set { collectedStepwiseOutputs }

    mergeResidualsOut = MergeResiduals(collectedStepwiseOutputs[0])
    prepareConditionalGenesisOut = PrepareConditionalGenesis(mergeResidualsOut[0])
    runConditionalGenesisOut = RunConditionalGenesis(prepareConditionalGenesisOut[0])

    runConditionalGenesisOut.collect().set { collectedConditionalOutputs }
    mergeConditionalGenesisOut = MergeConditionalGenesis(collectedConditionalOutputs[0])
}

