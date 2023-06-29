module Command

using JSON

using FASTQ

function concatenate_fastq(dna_read_directory, read_name_scheme)

    re_ = FASTQ.Raw.find(dna_read_directory)

    FASTQ.Raw.concatenate(re_, read_name_scheme)

end

function call_variants_on_bulk_cdna(
    output_directory,
    cdna_read_directory,
    number_of_jobs,
    reference_genome,
    molecule,
    exome,
    chromosome_position,
    chromosome_name,
    memory,
    tool_directory,
    snpeff,
    annotate_with_rsid,
    variant_database,
)

    for fi in (reference_genome, chromosome_position, chromosome_name, snpeff, variant_database)

        BioLab.Path.error_missing(fi)

    end

    pa = joinpath(output_directory, "call_variants_on_bulk_cdna")

    FASTQ.Support.error_if_directory(pa)

    pac = joinpath(pa, "align_cdna/")

    FASTQ.Support.error_if_directory(pac)

    re_ = FASTQ.Raw.find(cdna_read_directory)

    FASTQ.Raw.check_read(joinpath(pa, "check_read"), re_, number_of_jobs)

    FASTQ.Raw.align_cdna(pac, cdna_read_directory, reference_genome, number_of_jobs, al = "genome")

    pav = joinpath(pa, "call_germline_variant")

    for (ro, di_, fi_) in walkdir(pac)

        for fi in fi_
            if endswith(fi, ".bam")

                ba = joinpath(ro, fi)

                sa = basename(splitdir(ba)[1])

                pas = joinpath(pav, sa)

                FASTQ.BAM.call_germline_variant(
                    pas,
                    tool_directory,
                    ba,
                    reference_genome,
                    chromosome_position,
                    exome,
                    molecule,
                    number_of_jobs,
                    memory,
                    chromosome_name,
                    snpeff,
                    annotate_with_rsid,
                    variant_database,
                )

            end

        end

    end

end


function measure_gene_expression_of_bulk_cdna(
    output_directory,
    cdna_read_directory,
    number_of_jobs,
    reference_transcriptome,
    fragment_length,
    fragment_length_standard_deviation,
    organism,
    mouse_transcript_to_mouse_gene,
)

    for fi in (reference_transcriptome, mouse_transcript_to_mouse_gene)

        BioLab.Path.error_missing(fi)

    end

    pa = joinpath(output_directory, "measure_gene_expression_of_bulk_cdna")

    FASTQ.Support.error_if_directory(pa)

    pap = joinpath(pa, "psuedoalign/")

    FASTQ.Support.error_if_directory(pap)

    re_ = FASTQ.Raw.find(cdna_read_directory)

    FASTQ.Raw.check_read(joinpath(pa, "check_read"), re_, number_of_jobs)

    FASTQ.Raw.align_cdna(
        pap,
        cdna_read_directory,
        reference_transcriptome,
        number_of_jobs,
        al = "transcriptome",
        fr = fragment_length,
        sd = fragment_length_standard_deviation,
    )

    FASTQ.Abundance.make_gene_by_sample(pap, pa, organism, mouse_transcript_to_mouse_gene)

end

function measure_gene_expression_of_single_cell_cdna()

    for fi in (reference_genome,)

        BioLab.Path.error_missing(fi)

    end

    pa = joinpath(output_directory, "measure_gene_expression_of_single_cell_cdna")

    FASTQ.Support.error_if_directory(pa)

end

function call_variants_on_germline_dna(
    output_directory,
    read1,
    read2,
    number_of_jobs,
    memory,
    sample,
    reference_genome,
    chromosome_position,
    chromosome_name,
    snpeff,
    molecule,
    exome,
    tool_directory,
    annotate_with_rsid,
    variant_database,
)

    pa = joinpath(output_directory, "call_variants_on_germline_dna")

    FASTQ.Support.error_if_directory(pa)

    for fi in (read1, read2, reference_genome, chromosome_position, chromosome_name, snpeff)

        BioLab.Path.error_missing(fi)

    end

    FASTQ.Raw.check_read(pa, read1, read2, number_of_jobs)

    tr = joinpath(pa, "trim/")

    FASTQ.Raw.trim(tr, read1, read2, number_of_jobs)

    r1t = joinpath(tr, FASTQ.TRIMMED_R1)

    r2t = joinpath(tr, FASTQ.TRIMMED_R2)

    FASTQ.Raw.check_read(joinpath(pa, "check_trim"), [r1t, r2t], number_of_jobs)

    al = joinpath(pa, "align_dna")

    ba = joinpath(al, "$sample.bam")

    FASTQ.Raw.align_dna(al, sample, ba, r1t, r2t, reference_genome, number_of_jobs, memory)

    pav = joinpath(pa, "call_germline_variant")

    FASTQ.BAM.call_germline_variant(
        pav,
        tool_directory,
        ba,
        reference_genome,
        chromosome_position,
        exome,
        molecule,
        number_of_jobs,
        memory,
        chromosome_name,
        snpeff,
        annotate_with_rsid,
        variant_database,
    )

end

function call_variants_on_somatic_dna(
    output_directory,
    read1,
    read2,
    somatic_read1,
    somatic_read2,
    number_of_jobs,
    memory,
    sample,
    reference_genome,
    chromosome_position,
    chromosome_name,
    snpeff,
    molecule,
    exome,
    tool_directory,
    annotate_with_rsid,
    variant_database,
)

    pa = joinpath(output_directory, "call_variants_on_somatic_dna")

    FASTQ.Support.error_if_directory(pa)

    for fi in (
        read1,
        read2,
        somatic_read1,
        somatic_read2,
        reference_genome,
        chromosome_position,
        chromosome_name,
        snpeff,
    )

        BioLab.Path.error_missing(fi)

    end

    FASTQ.Raw.check_read(pa, read1, read2, number_of_jobs, somatic_read1, somatic_read2)

    trge = joinpath(pa, "trim", "germline")

    gr1 = joinpath(trge, FASTQ.TRIMMED_R1)

    gr2 = joinpath(trge, FASTQ.TRIMMED_R2)

    trso = joinpath(pa, "trim", "somatic")

    sr1 = joinpath(trso, FASTQ.TRIMMED_R1)

    sr2 = joinpath(trso, FASTQ.TRIMMED_R2)

    for g in [[trge, read1, read2], [trso, somatic_read1, somatic_read2]]

        FASTQ.Raw.trim(g[1], g[2], g[3], number_of_jobs)

    end

    FASTQ.Raw.check_read(joinpath(pa, "check_trim"), [gr1, gr2, sr1, sr2], number_of_jobs)

    alg = joinpath(pa, "align_$(molecule)_germline")

    als = joinpath(pa, "align_$(molecule)_somatic")

    bage = joinpath(alg, "$sample.bam")

    baso = joinpath(als, "$sample.bam")

    for g in [[alg, bage, gr1, gr2], [als, baso, sr1, sr2]]

        FASTQ.Raw.align_dna(
            g[1],
            sample,
            g[2],
            g[3],
            g[4],
            reference_genome,
            number_of_jobs,
            memory,
        )

    end

    pav = joinpath(pa, "call_somatic_variant")

    bagem = joinpath(alg, "$sample.bam")

    basom = joinpath(als, "$sample.bam")

    FASTQ.BAM.call_somatic_variant(
        pav,
        tool_directory,
        bagem,
        reference_genome,
        chromosome_position,
        basom,
        exome,
        number_of_jobs,
        memory,
        chromosome_name,
        snpeff,
        annotate_with_rsid,
        variant_database,
    )

end

function benchmark(
    output_directory,
    reference_genome,
    rtg_tools,
    number_of_jobs,
    name_chromosome,
    query_vcf,
    truth_vcf,
    confident_regions_bed,
)

    FASTQ.Support.log()

    pa = joinpath(output_directory, "benchmark")

    FASTQ.Support.error_if_directory(pa)

    for fi in (name_chromosome, query_vcf, truth_vcf, confident_regions_bed)

        BioLab.Path.error_missing(fi)

    end

    red = split(reference_genome, ".gz")[1]

    if !isfile(red)

        run(`bgzip -d $reference_genome`)

    end

    sd = replace(red, "fna" => "sdf")

    if !isdir(sd)

        @info "Making vcfeval genome sdf"

        run(`$rtg_tools format -o $sd $red`)

    end

    vqn = replace(query_vcf, "pass" => "pass_rename_chromosomes")

    if !isfile(vqn)

        @warn "Renaming query VCF chromosomes"

        run(
            `bcftools annotate --threads=$number_of_jobs --rename-chrs=$name_chromosome --output=$vqn $query_vcf`,
        )

    end

    @info "Running vcfeval"

    ouv = joinpath(pa, "vcfeval")

    rte = joinpath(rtg_tools, "rtg")

    run(`$rte vcfeval 
        --baseline=$truth_vcf 
        --bed-regions=$confident_regions_bed 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$number_of_jobs`)

    @info "Running hap.py"

    ouh = joinpath(pa, "happy/")

    mkdir(ouh)

    ho = "/home/"

    vouh = joinpath(ho, splitpath(ouh)[end])

    pvt = dirname(BioLab.Path.make_absolute(truth_vcf))

    vvt = joinpath(ho, basename(pvt))

    pvqn = dirname(BioLab.Path.make_absolute(vqn))

    vvqn = joinpath(ho, basename(pvqn))

    pbd = dirname(BioLab.Path.make_absolute(confident_regions_bed))

    vbd = joinpath(ho, "confident_regions_bed/")

    pre = dirname(BioLab.Path.make_absolute(red))

    vre = joinpath(ho, basename(pre))

    vrt = joinpath(ho, basename(rtg_tools))

    vsd = joinpath(ho, basename(sd))

    id = readlines(pipeline(`docker run 
            --interactive 
            --detach 
            --tty 
            --user root
            -v $pvt:$vvt 
            -v $pvqn:$vvqn 
            -v $pbd:$vbd 
            -v $pre:$vre 
            -v $ouh:$vouh 
            -v $rtg_tools:$vrt 
            -v $sd:$vsd 
            pkrusche/hap.py
            bash`))


    vh = joinpath(vouh, "hap.py")

    vr = joinpath(vre, basename(red))

    vb = joinpath(vbd, basename(confident_regions_bed))

    vtr = joinpath(vvt, basename(truth_vcf))

    vqu = joinpath(vvqn, basename(vqn))

    readlines(
        pipeline(
            `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $vtr $vqu -f $vb -r $vr -o $vh --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
        ),
    )

    FASTQ.Support.remove_docker_container(id)

end

end
