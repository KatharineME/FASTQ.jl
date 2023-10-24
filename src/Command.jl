module Command

using ..FASTQ

const RN = "R1"

function concatenate_fastq(dna_read_directory; read_name_scheme = RN)

    FASTQ.Support.log_top_level_function()

    di_ = readdir(dna_read_directory, join = true)

    for sa in readdir(dna_read_directory, join = true)

        if isdir(sa)

            fq_ = FASTQ.Raw.find(sa)

            if length(fq_) > 2

                pa = FASTQ.Support.trash_remake_directory(
                    joinpath(splitdir(sa)[1], string(basename(sa), "Concatenated")),
                )

                FASTQ.Raw.concatenate(pa, fq_; na = read_name_scheme)

            end

        end

    end

    nothing

end

function _get_snpeff_path(to)

    joinpath(to, "snpEff", "snpEff.jar")

end

function _combine_and_annotate_vcf(pase, pass, paco, vc_, re, chn, va, se, n_jo, me; sa = nothing)

    if sa != nothing

        pase, pass = joinpath(pase, sa), joinpath(pass, sa)

    end

    FASTQ.VCF.combine_vcf(paco, vc_, chn, n_jo)

    papa = FASTQ.VCF.annotate_with_snpeff(pase, paco, re, se, n_jo, me)

    FASTQ.VCF.annotate_with_snpsift(pass, papa, va, se, n_jo)

end

const CTR = "CheckTrim"

const TR = "Trim"

const ALD = "AlignDNA"

const CG = "CallGermlineVariant"

const AN = "Annotate"

const SE = "Snpeff"

const SS = "Snpsift"

const CRA = "1.CheckRaw"

const TW = "2."

const TH = "3."

const FO = "4."

const FI = "5."

const SI = "6."

function call_variants_on_germline_dna(
    output_directory,
    dna_read_directory,
    exome,
    reference_genome,
    variant_database,
    tool_directory,
    number_of_jobs,
    memory;
    read_name_scheme = RN,
)

    FASTQ.Support.log_top_level_function()

    chs, chn = FASTQ.Reference.get_chromosome_file_path(reference_genome)

    se = _get_snpeff_path(tool_directory)

    FASTQ.Support.error_if_file_missing((reference_genome, chs, chn, se, variant_database))

    sa_fq_ = FASTQ.Support.make_sample_to_fastq_dictionary(dna_read_directory, read_name_scheme)

    cha, tr, cht, al, va, pase, pass = FASTQ.Support.make_analysis_directory(
        output_directory,
        "CallVariantsonGermlineDNA",
        (
            CRA,
            string(TW, TR),
            string(TH, CTR),
            string(FO, ALD),
            string(FI, CG),
            joinpath(string(SI, AN), SE),
            joinpath(string(SI, AN), SS),
        );
        sa_fq_ = sa_fq_,
    )

    FASTQ.Reference.index_genome_file(reference_genome, chs)

    for (sa, fq_) in sa_fq_

        san, fq1, fq2 = FASTQ.Raw.check(joinpath(cha, basename(sa)), fq_, number_of_jobs)

        r1t, r2t = FASTQ.Raw.trim(joinpath(tr, san), fq1, fq2, number_of_jobs)

        FASTQ.Raw.check(joinpath(cht, san), (r1t, r2t), number_of_jobs)

        ba = FASTQ.Raw.align_dna(
            joinpath(al, san),
            san,
            r1t,
            r2t,
            reference_genome,
            number_of_jobs,
            memory,
        )

        vc_, paco = FASTQ.BAM.call_germline_variant(
            joinpath(va, san),
            ba,
            "dna",
            exome,
            reference_genome,
            chs,
            chn,
            variant_database,
            tool_directory,
            number_of_jobs,
            memory,
        )

        _combine_and_annotate_vcf(
            pase,
            pass,
            paco,
            vc_,
            reference_genome,
            chn,
            variant_database,
            se,
            number_of_jobs,
            memory;
            sa = san,
        )

    end

    nothing

end

const GE = "Germline"

const SO = "Somatic"

function call_variants_on_somatic_dna(
    output_directory,
    read1,
    read2,
    somatic_read1,
    somatic_read2,
    exome,
    sample,
    reference_genome,
    variant_database,
    tool_directory,
    number_of_jobs,
    memory,
)

    FASTQ.Support.log_top_level_function()

    chs, chn = FASTQ.Reference.get_chromosome_file_path(reference_genome)

    se = _get_snpeff_path(tool_directory)

    FASTQ.Support.error_if_file_missing((
        read1,
        read2,
        somatic_read1,
        somatic_read2,
        reference_genome,
        chs,
        chn,
        se,
    ))

    tr, ald, an = string(TW, TR), string(FO, ALD), string(SI, AN)

    cha, trge, trso, cht, alg, als, pav, pase, pass = FASTQ.Support.make_analysis_directory(
        output_directory,
        "CallVariantsonSomaticDNA",
        (
            CRA,
            joinpath(tr, GE),
            joinpath(tr, SO),
            string(TH, CTR),
            joinpath(ald, GE),
            joinpath(ald, SO),
            "5.CallSomaticVariant",
            joinpath(an, SE),
            joinpath(an, SS),
        ),
    )

    FASTQ.Raw.check(cha, (read1, read2, somatic_read1, somatic_read2), number_of_jobs)

    gr1, gr2 = FASTQ.Raw.trim(trge, read1, read2, number_of_jobs)

    sr1, sr2 = FASTQ.Raw.trim(trso, somatic_read1, somatic_read2, number_of_jobs)

    FASTQ.Raw.check(cht, (gr1, gr2, sr1, sr2), number_of_jobs)

    ba_ = Vector{String}()

    for (pa, r1, r2) in ((alg, gr1, gr2), (als, sr1, sr2))

        push!(
            ba_,
            FASTQ.Raw.align_dna(pa, sample, r1, r2, reference_genome, number_of_jobs, memory),
        )

    end

    bagem, basom = ba_[1], ba_[2]

    FASTQ.Reference.index_genome_file(reference_genome, chs)

    vc_, paco = FASTQ.BAM.call_somatic_variant(
        pav,
        bagem,
        basom,
        exome,
        reference_genome,
        chs,
        chn,
        variant_database,
        tool_directory,
        number_of_jobs,
        memory,
    )

    _combine_and_annotate_vcf(
        pase,
        pass,
        paco,
        vc_,
        reference_genome,
        chn,
        variant_database,
        se,
        number_of_jobs,
        memory,
    )

    nothing

end

function call_variants_on_bulk_cdna(
    output_directory,
    cdna_read_directory,
    exome,
    reference_genome,
    variant_database,
    tool_directory,
    number_of_jobs,
    memory;
    read_name_scheme = RN,
    gene_annotation = nothing,
)

    FASTQ.Support.log_top_level_function()

    chs, chn = FASTQ.Reference.get_chromosome_file_path(reference_genome)

    se = _get_snpeff_path(tool_directory)

    FASTQ.Support.error_if_file_missing((reference_genome, chs, chn, se, variant_database))

    sa_fq_ = FASTQ.Support.make_sample_to_fastq_dictionary(cdna_read_directory, read_name_scheme)

    cha, al, va, pase, pass = FASTQ.Support.make_analysis_directory(
        output_directory,
        "CallVariantsonBulkCDNA",
        (
            CRA,
            string(TW, "AlignBulkCDNAtoGenome"),
            string(TH, CG),
            joinpath(string(FO, AN), SE),
            joinpath(string(FO, AN), SS),
        );
        sa_fq_ = sa_fq_,
    )

    id = FASTQ.Reference.generate_star_genome_file(
        reference_genome,
        number_of_jobs;
        ga = gene_annotation,
    )

    for (sa, fq_) in sa_fq_

        san, fq1, fq2 = FASTQ.Raw.check(joinpath(cha, basename(sa)), fq_, number_of_jobs)

        ba = FASTQ.Raw.align_bulk_cdna_to_genome(joinpath(al, san), fq1, fq2, id, number_of_jobs)

        vc_, paco = FASTQ.BAM.call_germline_variant(
            joinpath(va, san),
            ba,
            "cdna",
            exome,
            reference_genome,
            chs,
            chn,
            variant_database,
            tool_directory,
            number_of_jobs,
            memory,
        )

        _combine_and_annotate_vcf(
            pase,
            pass,
            paco,
            vc_,
            reference_genome,
            chn,
            variant_database,
            se,
            number_of_jobs,
            memory;
            sa = san,
        )

    end

    nothing

end

function measure_gene_expression_of_bulk_cdna(
    output_directory,
    cdna_read_directory,
    organism,
    fragment_length,
    fragment_length_standard_deviation,
    reference,
    number_of_jobs;
    read_name_scheme = RN,
    method = "align_to_transcriptome",
    gene_annotation = nothing,
    mouse_transcript_to_mouse_gene = nothing,
)

    FASTQ.Support.log_top_level_function()

    FASTQ.Support.error_if_file_missing((reference,))

    sa_fq_ = FASTQ.Support.make_sample_to_fastq_dictionary(cdna_read_directory, read_name_scheme)

    me_na_ = Dict(
        "align_to_transcriptome" => "2.AlignBulkCDNAtoTranscriptome",
        "align_to_genome" => "2.AlignandQuantifyBulkCDNAtoGenome",
    )

    cha, al = FASTQ.Support.make_analysis_directory(
        output_directory,
        "MeasureGeneExpressionofBulkCDNA",
        (CRA, me_na_[method]);
        sa_fq_ = sa_fq_,
    )

    if method == "align_to_transcriptome"

        for (sa, fq_) in sa_fq_

            san, fq1, fq2 = FASTQ.Raw.check(joinpath(cha, basename(sa)), fq_, number_of_jobs)

            FASTQ.Raw.align_bulk_cdna_to_transcriptome(
                joinpath(al, san),
                fq1,
                fq2,
                fragment_length,
                fragment_length_standard_deviation,
                reference,
                number_of_jobs,
            )

            #FASTQ.Abundance.make_gene_by_sample(pal, wo, organism, mouse_transcript_to_mouse_gene)

        end

    elseif method == "align_to_genome"

        id = FASTQ.Reference.generate_star_genome_file(
            reference,
            number_of_jobs;
            ga = gene_annotation,
        )

        for (sa, fq_) in sa_fq_

            san, fq1, fq2 = FASTQ.Raw.check(joinpath(cha, basename(sa)), fq_, number_of_jobs)

            ba = FASTQ.Raw.align_and_quantify_bulk_cdna_to_genome(
                joinpath(al, san),
                fq1,
                fq2,
                id,
                number_of_jobs,
            )

            #FASTQ.Abundance.make_gene_by_sample()

        end

    end

    nothing

end

function measure_gene_expression_of_single_cell_cdna(
    output_directory,
    single_cell_read_directory,
    reference_genome,
    number_of_jobs;
    read_name_scheme = RN,
    gene_annotation = nothing,
    whitelist = nothing,
    barcodestart = 1,
    barcodelength = 16,
    umistart = 17,
    umilength = 12,
    readlength = 151,
)

    FASTQ.Support.log_top_level_function()

    FASTQ.Support.error_if_file_missing((reference_genome,))

    sa_fq_ =
        FASTQ.Support.make_sample_to_fastq_dictionary(single_cell_read_directory, read_name_scheme)

    id = FASTQ.Reference.generate_star_genome_file(
        reference_genome,
        number_of_jobs;
        ga = gene_annotation,
    )

    cha, al = FASTQ.Support.make_analysis_directory(
        output_directory,
        "MeasureGeneExpressionofSingleCellCDNA",
        (CRA, "2.AlignSingleCellCDNAtoGenome");
        sa_fq_ = sa_fq_,
    )

    for (sa, fq_) in sa_fq_

        san, fq1, fq2 = FASTQ.Raw.check(joinpath(cha, basename(sa)), fq_, number_of_jobs)

        FASTQ.Raw.align_single_cell_cdna_to_genome(
            joinpath(al, san),
            fq1,
            fq2,
            id,
            number_of_jobs;
            wh = whitelist,
            bas = barcodestart,
            bal = barcodelength,
            ums = umistart,
            uml = umilength,
            rel = readlength,
        )

    end

    nothing

end

end
