function check_read(fq_, di, n_jo)

    @assert make_directory(di, "check read")

    println("Running FastQC ...")

    run(`fastqc --threads $(minimum((length(fq_), n_jo))) --outdir $di $fq_`)

    println("Running MultiQC ...")

    run(`multiqc --outdir $di $di`)

    return

end
