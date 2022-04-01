function check_read(re_, di, n_jo)

    Fastq.support.log()

    Fastq.support.error_if_directory(di)

    run(`fastqc --threads $(minimum((length(re_), n_jo))) --outdir $di $re_`)

    run(`multiqc --outdir $di $di`)

end
