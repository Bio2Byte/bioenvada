process compressDirectory {
    publishDir "$params.outDir", mode: "copy"

    input:
    val compressedFileName
    path outputFiles

    output:
    path "${compressedFileName}.tar.gz"
    val true ,emit:gate

    script:
    """
    echo $compressedFileName
    tar -czh -f ${compressedFileName}.tar.gz $outputFiles
    """
}