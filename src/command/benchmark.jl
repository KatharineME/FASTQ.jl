function benchmark(se)

    fe_va = read_setting(se)

    pou, r1, r2, n_jo, me, sa, ge, chs, chn, sn = fe_va["output_directory"],
    fe_va["truth_vcf"],
    fe_va["query_vcf"],
    fe_va["confident_regions_bed"],
    fe_va["rtg_tools"],
    fe_va["reference_genome"],
    fe_va["chromosome_position"],
    fe_va["chromosome_name"],
    fe_va["snpeff"]

    pa = joinpath(pou, "benchmark")

    Fastq.support.error_if_directory(pa)



end
