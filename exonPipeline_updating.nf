params.trait = 'BMI'
params.numSlice = 1
params.minBinWidth = 45
params.minCPM = 3
params.minPercentsExpressed = 0.015
params.number_of_plates = 28
params.visitcode = "v1"
params.binLevelPval = 0.01
params.binsRDSpath = '/project/data/bins.RDS'
params.hgnc2ensemblMapPath = '/project/data/HGNC2ensembl.RDS'

nextflow.enable.dsl=2

process PreprocessData {
    container 'edkang0925/exonproject-m1:latest'

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
    container 'edkang0925/exonproject-m1:latest'

    input:
    path(pc)
    path(vst)

    output:
    path("outputs/beforeRegression/${params.trait}/sliced/*")
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
    container 'edkang0925/exonproject-m1:latest'

    input:
    path(slicedBinCount)
    path(phenotypeBeforeRegression)
    path(combinedDF)

    output:
    path("outputs/afterRegression/${params.trait}/sliced/*")

    script:
    """
    Rscript /project/scripts/stepwise_regression_perbin.R \
        --combined_df ${combinedDF} \
        --exon_expression_path ${slicedBinCount} \
        --output_directory "outputs/afterRegression/${params.trait}/sliced/"
    """
}

process MergeResiduals {
    container 'edkang0925/exonproject-m1:latest'
    publishDir "./", mode: 'copy'

    input:
    path x

    output:
    path("outputs/afterRegression/${params.trait}/combined/residualExonExpression.RDS")

    script:
    """
    Rscript /project/scripts/combine_residuals.R \
        --input_files "${x.join(',')}" \
        --output_prefix "outputs/afterRegression/${params.trait}/combined/"
    """
}

process PrepareConditionalGenesis {
    container 'edkang0925/exonproject-m1:latest'

    input:
    path(mergedResidualsExpression)

    output:
    path("outputs/conditionalGenesisInput/${params.trait}/conditionalSlice_*.RDS")

    script:
    """
    Rscript /project/scripts/preprocessForConditional_all_updating.R \
        --trait ${params.trait} \
        --numSlice ${params.numSlice} \
        --bins ${params.binsRDSpath} \
        --hgnc2ensembl ${params.hgnc2ensemblMapPath} \
        --mergedResidualsPath ${mergedResidualsExpression} \
        --outputDir "outputs/conditionalGenesisInput/${params.trait}/"
    """
}

process RunConditionalGenesis {
    container 'edkang0925/exonproject-m1:latest'

    input:
    path(conditionalGenesisInputSliceFile)

    output:
    path("outputs/conditionalGenesisOutput/${params.trait}/conditionalGenesisResultSlice_*.RDS")

    script:
    """
    Rscript /project/scripts/conditionalGenesis_chunk.R \
        --trait ${params.trait} \
        --exon.path ${conditionalGenesisInputSliceFile} \
        --outputDir "outputs/conditionalGenesisOutput/${params.trait}/"
    """
}

process MergeConditionalGenesis {
    container 'edkang0925/exonproject-m1:latest'
    publishDir "./", mode: 'copy'

    input:
    path x

    output:
    path("outputs/combinedConditionalGenesis/${params.trait}/comnbinedConditionalGenesis_${params.trait}.RDS")

    script:
    """
    Rscript /project/scripts/combineConditionalGenesis.R \
        --trait ${params.trait} \
        --input_files ${x.join(',')} \
        --outputFile "outputs/combinedConditionalGenesis/${params.trait}/comnbinedConditionalGenesis_${params.trait}.RDS"
    """
}


workflow {
    preprocessDataOut = PreprocessData()

    prepareRegressionOut = PrepareRegression(preprocessDataOut[0], preprocessDataOut[1])

    stepwiseRegressionOut =StepwiseRegression(prepareRegressionOut[0]|flatten, prepareRegressionOut[1], prepareRegressionOut[2])
    // Collect StepwiseRegressionOut as MergeResiduals need all of the outputs in the dir

    mergeResidualsOut = MergeResiduals(stepwiseRegressionOut.collect())
    prepareConditionalGenesisOut = PrepareConditionalGenesis(mergeResidualsOut[0])
    runConditionalGenesisOut = RunConditionalGenesis(prepareConditionalGenesisOut[0]|flatten)
    // Collect RunConditionalGenesisOut before merging them
    mergeConditionalGenesisOut = MergeConditionalGenesis(runConditionalGenesisOut.collect())
}

