params.trait = 'BMI'
params.numSlice = 1
params.minBinWidth = 45
params.minCPM = 3
params.minPercentsExpressed = 0.015
params.number_of_plates = 28
params.visitcode = "v1"
params.binLevelPval = 0.01

nextflow.enable.dsl=2

process PreprocessData {
    container 'edkang0925/exonproject-m1:latest'

    output:
    val(params.trait)
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
    container 'edkang0925/exonproject-m1:latest'

    input:
    val(trait)
    path(pc)
    path(vst)

    output:
    val(trait)
    path("outputs/beforeRegression/${trait}/sliced/*")
    path("outputs/beforeRegression/${trait}/phenotype_df_beforeRegression_${trait}.RDS")
    path("outputs/beforeRegression/${trait}/combined_df_beforeRegression_${trait}.RDS")

    script:
    """
    Rscript /project/scripts/prepareRegression.R --trait ${trait} --numSlice ${params.numSlice} \
        --number_of_plates ${params.number_of_plates} --visitcode ${params.visitcode} \
        --pcFilePath ${pc} --vstFilePath ${vst}
    """
}

process StepwiseRegression {
    container 'edkang0925/exonproject-m1:latest'

    input:
    val(trait)
    path(slicedBinCount)
    path(phenotypeBeforeRegression)
    path(combinedDF)

    output:
    val(trait)
    path("outputs/afterRegression/${trait}/sliced/")

    script:
    """
    Rscript /project/scripts/stepwise_regression_perbin.R \
        --combined_df ${combinedDF} \
        --exon_expression_path ${slicedBinCount} \
        --output_directory "outputs/afterRegression/${trait}/sliced/"
    """
}

process MergeResiduals {
    container 'edkang0925/exonproject-m1:latest'
    publishDir "./", mode: 'copy'

    input:
    val(trait)
    path(residualBinExpression)

    output:
    val(trait)
    path("outputs/afterRegression/${trait}/combined/residualExonExpression.RDS")

    script:
    """
    Rscript /project/scripts/combine_residuals.R \
        --input_prefix ${residualBinExpression} \
        --output_prefix "outputs/afterRegression/${trait}/combined/"
    """
}

process PrepareConditionalGenesis {
    container 'edkang0925/exonproject-m1:latest'

    input:
    val(trait)
    path(mergedResidualsExpression)

    output:
    val(trait)
    path("outputs/conditionalGenesisInput/${trait}/*")

    script:
    """
    Rscript /project/scripts/preprocessForConditional_all.R \
        --trait ${trait} \
        --numSlice ${params.numSlice} \
        --pvalThreshold ${params.binLevelPval} \
        --mergedResidualsPath ${mergedResidualsExpression} \
        --outputDir "outputs/conditionalGenesisInput/${trait}/"
    """
}

process RunConditionalGenesis {
    container 'edkang0925/exonproject-m1:latest'

    input:
    val(trait)
    path(conditionalGenesisInputSliceFile)

    output:
    val(trait)
    path("outputs/conditionalGenesisOutput/${trait}/")

    script:
    """
    Rscript /project/scripts/conditionalGenesis_chunk.R \
        --trait ${trait} \
        --exon.path ${conditionalGenesisInputSliceFile} \
        --outputFile "outputs/conditionalGenesisOutput/${trait}/results_001.RDS"
    """
}

process MergeConditionalGenesis {
    container 'edkang0925/exonproject-m1:latest'
    publishDir "./", mode: 'copy'

    input:
    val(trait)
    path(conditionalGenesisOutputDir)

    output:
    val(trait)
    path("outputs/combinedConditionalGenesis/${trait}/comnbinedConditionalGenesis_${trait}.RDS")

    script:
    """
    Rscript /project/scripts/combineConditionalGenesis.R \
        --trait ${trait} \
        --inputDir ${conditionalGenesisOutputDir} \
        --outputFile "outputs/combinedConditionalGenesis/${trait}/comnbinedConditionalGenesis_${trait}.RDS"
    """
}


workflow {
    preprocessDataOut = PreprocessData()
    prepareRegressionOut = PrepareRegression(preprocessDataOut[0], preprocessDataOut[1], preprocessDataOut[2])
    stepwiseRegressionOut = StepwiseRegression(prepareRegressionOut[0], prepareRegressionOut[1], prepareRegressionOut[2], prepareRegressionOut[3])
    mergeResidualsOut = MergeResiduals(stepwiseRegressionOut[0], stepwiseRegressionOut[1])
    prepareConditionalGenesisOut = PrepareConditionalGenesis(mergeResidualsOut[0], mergeResidualsOut[1])
    runConditionalGenesisOut = RunConditionalGenesis(prepareConditionalGenesisOut[0], prepareConditionalGenesisOut[1])
    mergeConditionalGenesisOut = MergeConditionalGenesis(runConditionalGenesisOut[0], runConditionalGenesisOut[1])
}

