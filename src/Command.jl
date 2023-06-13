module Command

using JSON

using FASTQ

function concatenate_fastq(dna_read_directory, read_name_scheme)

    re_ = Fastq.fastq.find(dna_read_directory)

    Fastq.fastq.concatenate(re_, read_name_scheme)

end

function apply_cdna_to_genome(
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

    pou = joinpath(output_directory, "apply_cdna_to_genome")

    Fastq.support.error_if_directory(pou)

    pac = joinpath(pou, "align_cdna/")

    Fastq.support.error_if_directory(pac)

    re_ = Fastq.fastq.find(cdna_read_directory)

    Fastq.fastq.check_read(re_, joinpath(pou, "check_read"), number_of_jobs)

    Fastq.fastq.align_cdna_samples(
        pac,
        cdna_read_directory,
        reference_genome,
        number_of_jobs,
        al = "genome",
    )

    pav = joinpath(pou, "call_germline_variant")

    for (ro, di_, fi_) in walkdir(pac)

        for fi in fi_
            if endswith(fi, ".bam")

                ba = joinpath(ro, fi)

                sa = basename(splitdir(ba)[1])

                pas = joinpath(pav, sa)

                Fastq.bam.call_germline_variant(
                    molecule,
                    exome,
                    ba,
                    reference_genome,
                    chromosome_position,
                    chromosome_name,
                    pas,
                    number_of_jobs,
                    memory,
                    tool_directory,
                    snpeff,
                    annotate_with_rsid,
                    variant_database,
                )

            end

        end

    end

end


function apply_cdna_to_transcriptome(
    output_directory,
    cdna_read_directory,
    number_of_jobs,
    reference_transcriptome,
    fragment_length,
    fragment_length_standard_deviation,
    organism,
    mouse_transcript_to_mouse_gene,
)

    pou = joinpath(output_directory, "apply_cdna_to_transcriptome")

    Fastq.support.error_if_directory(pou)

    pap = joinpath(pou, "psuedoalign/")

    Fastq.support.error_if_directory(pap)

    re_ = Fastq.fastq.find(cdna_read_directory)

    Fastq.fastq.check_read(re_, joinpath(pou, "check_read"), number_of_jobs)

    Fastq.fastq.align_cdna_samples(
        pap,
        cdna_read_directory,
        reference_transcriptome,
        number_of_jobs,
        al = "transcriptome",
        fr = fragment_length,
        sd = fragment_length_standard_deviation,
    )

    Fastq.abundance.make_gene_by_sample(pap, pou, organism, mouse_transcript_to_mouse_gene)

end

function apply_germline_dna_to_genome(
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

    pa = joinpath(output_directory, "apply_germline_dna_to_genome")

    Fastq.support.error_if_directory(pa)

    Fastq.fastq.examine_read(read1, read2, pa, number_of_jobs)

    for pa in [read1, read2, reference_genome, chromosome_position, chromosome_name, snpeff]

        if !isfile(pa)

            error("$pa does not exist.")

        end

    end

    tr = joinpath(pa, "trim/")

    Fastq.fastq.trim(read1, read2, tr, number_of_jobs)

    r1t = joinpath(tr, Fastq.TRIMMED_R1)

    r2t = joinpath(tr, Fastq.TRIMMED_R2)

    Fastq.fastq.check_read([r1t, r2t], joinpath(pa, "check_trim"), number_of_jobs)

    al = joinpath(pa, "align_dna")

    ba = joinpath(al, "$sample.bam")

    Fastq.fastq.align_dna(al, sample, ba, r1t, r2t, reference_genome, number_of_jobs, memory)

    pav = joinpath(pa, "call_germline_variant")

    Fastq.bam.call_germline_variant(
        molecule,
        exome,
        ba,
        reference_genome,
        chromosome_position,
        chromosome_name,
        pav,
        number_of_jobs,
        memory,
        tool_directory,
        snpeff,
        annotate_with_rsid,
        variant_database,
    )

end

function apply_somatic_dna_to_genome(
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

    pa = joinpath(output_directory, "apply_somatic_dna_to_genome")

    Fastq.support.error_if_directory(pa)

    Fastq.fastq.examine_read(read1, read2, pa, number_of_jobs, somatic_read1, somatic_read2)

    for fi in [
        read1,
        read2,
        somatic_read1,
        somatic_read2,
        reference_genome,
        chromosome_position,
        chromosome_name,
        snpeff,
    ]

        if !isfile(fi)

            error("$fi does not exist.")

        end

    end

    trge = joinpath(pa, "trim", "germline")

    gr1 = joinpath(trge, Fastq.TRIMMED_R1)

    gr2 = joinpath(trge, Fastq.TRIMMED_R2)

    trso = joinpath(pa, "trim", "somatic")

    sr1 = joinpath(trso, Fastq.TRIMMED_R1)

    sr2 = joinpath(trso, Fastq.TRIMMED_R2)

    for g in [[read1, read2, trge], [somatic_read1, somatic_read2, trso]]

        Fastq.fastq.trim(g[1], g[2], g[3], number_of_jobs)

    end

    Fastq.fastq.check_read([gr1, gr2, sr1, sr2], joinpath(pa, "check_trim"), number_of_jobs)

    alg = joinpath(pa, "align_$(molecule)_germline")

    als = joinpath(pa, "align_$(molecule)_somatic")

    bage = joinpath(alg, "$sample.bam")

    baso = joinpath(als, "$sample.bam")

    for g in [[alg, bage, gr1, gr2], [als, baso, sr1, sr2]]

        Fastq.fastq.align_dna(
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

    Fastq.bam.call_somatic_variant(
        exome,
        bagem,
        basom,
        reference_genome,
        chromosome_position,
        chromosome_name,
        pav,
        number_of_jobs,
        memory,
        tool_directory,
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

    Fastq.support.log()

    pa = joinpath(output_directory, "benchmark")

    Fastq.support.error_if_directory(pa)


    # Make vcfeval sdf

    red = split(reference_genome, ".gz")[1]

    if !isfile(reference_genome)

        run(`bgzip -d $reference_genome`)

    end

    sd = replace(red, "fna" => "sdf")

    if !isdir(sd)

        println("\nMaking vcfeval genome sdf\n")

        run(`$rtg_tools format -o $sd $red`)

    end


    # Rename chromosomes from numbers to strings

    vqn = replace(query_vcf, "pass" => "pass_rename_chromosomes")

    if !isfile(vqn)

        println("\nRenaming query VCF chromosomes\n")

        run(
            `bcftools annotate --threads=$number_of_jobs --rename-chrs=$name_chromosome --output=$vqn $query_vcf`,
        )

    end


    # Run vcfeval

    println("\nRunning vcfeval\n")

    ouv = joinpath(pa, "vcfeval")

    rte = joinpath(rtg_tools, "rtg")

    run(`$rte vcfeval 
        --baseline=$truth_vcf 
        --bed-regions=$confident_regions_bed 
        --calls=$vqn 
        --template=$sd 
        --output=$ouv
        --threads=$number_of_jobs`)


    # Run hap.py in container

    println("\nRunning hap.py\n")

    ouh = joinpath(pa, "happy/")

    mkdir(ouh)

    ho = "/home/"

    vouh = joinpath(ho, splitpath(ouh)[end])

    pvt = dirname(Fastq.support.get_full_path(truth_vcf))

    vvt = joinpath(ho, basename(pvt))

    pvqn = dirname(Fastq.support.get_full_path(vqn))

    vvqn = joinpath(ho, basename(pvqn))

    pbd = dirname(Fastq.support.get_full_path(confident_regions_bed))

    vbd = joinpath(ho, "confident_regions_bed/")

    pre = dirname(Fastq.support.get_full_path(red))

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

    readlines(
        pipeline(
            `docker exec --interactive $id bash -c "/opt/hap.py/bin/hap.py $(joinpath(vvt, basename(truth_vcf))) $(joinpath(vvqn, basename(vqn))) -f $(joinpath(vbd, basename(confident_regions_bed))) -r $(joinpath(vre, basename(red))) -o $(joinpath(vouh, "hap.py")) --engine-vcfeval-path $vrt --engine-vcfeval-template $vsd"`,
        ),
    )

    Fastq.support.remove_docker_container(id)

end

end
