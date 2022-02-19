function check_read(re_, di, n_jo)

    @assert make_directory(di, "check read")

    println("Running FastQC ...")

    run(`fastqc --threads $(minimum((length(re_), n_jo))) --outdir $di $re_`)

    println("Running MultiQC ...")

    run(`multiqc --outdir $di $di`)

    return

end
