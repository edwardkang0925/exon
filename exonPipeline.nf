params.trait = 'BMI'

process SplitData {
    container 'your-docker-image'

    output:
    file('split_data/${params.trait}/*') into splitDataCh

    script:
    """
    Rscript script1.R
    """
}

process ProcessSplitData {
    container 'your-docker-image'

    input:
    file(inputData) from splitDataCh

    output:
    file("${inputData}.out") into processedDataCh

    script:
    """
    Rscript script2.R ${inputData}
    """
}

process MergeData {
    container 'your-docker-image'

    input:
    file(inputData) from processedDataCh.collect()

    script:
    """
    Rscript script3.R
    """
}
