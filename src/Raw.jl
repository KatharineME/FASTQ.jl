module Raw

using FASTQ

function find(di)

    FASTQ.Support.log()

    re_ = []

    na_n = Dict(".fq" => 0, ".fastq" => 0, "fq.gz" => 0, "fastq.gz" => 0)

    for (ro, di_, fi_) in walkdir(di)

        for fi in fi_

            if !occursin(".md5", fi)

                for (na, _) in na_n

                    if endswith(fi, na)

                        na_n[na] += 1

                        if endswith(fi, ".gz")

                            push!(re_, joinpath(ro, fi))

                        end

                    end

                end

            end

        end

    end

    println("\nFile types found in $di:\n")

    for na in na_n
        println(na)
    end

    println("\nSize of gzipped files:\n")

    for fi in re_
        println("File $fi is: $(Base.format_bytes(stat(fi).size))")
    end

    return re_

end

function check_read(re_, di, n_jo)

    FASTQ.Support.log()

    FASTQ.Support.error_if_directory(di)

    run(`fastqc --threads $(minimum((length(re_), n_jo))) --outdir $di $re_`)

    run(`multiqc --outdir $di $di`)

end

function examine_read(r1, r2, pa, n_jo, sor1 = nothing, sor2 = nothing)

    if sor1 === nothing

        re_ = [r1, r2]

    else

        re_ = [r1, r2, sor1, sor2]

    end

    FASTQ.Raw.check_read(re_, joinpath(pa, "check_raw"), n_jo)

end

function concatenate(fq_, na = "R1")

    FASTQ.Support.log()

    fo_ = []

    re_ = []

    for fq in fq_

        if occursin(na, fq)

            push!(fo_, fq)

        elseif occursin(replace(na, "1" => "2"), fq)

            push!(re_, fq)

        end

    end

    n_fo = length(fo_)

    n_re = length(re_)

    println("Number of forward read files found = $n_fo\n")

    println("Number of reverse read files found = $n_re\n")

    sa = last(splitdir(dirname(fq_[1])))

    co = joinpath(dirname(dirname(fq_[1])), "$(sa)_concatenated")

    if n_fo <= 1 && n_re <= 1

        println(
            "\nNothing to concatenate, number of forward reads and reverse reads are both <= 1\n",
        )

    else

        FASTQ.Support.error_if_directory(co)

        println("\nConcatenating\n")

        gr_su = Dict(fo_ => "_R1.fastq.gz", re_ => "_R2.fastq.gz")

        for gr in keys(gr_su)

            run(pipeline(`cat $gr`, stdout = joinpath(co, "$(sa)$(gr_su[gr])")))

        end

        println("\nConcatenated files saved at: $co\n")

    end

end

function trim(r1, r2, pa, n_jo)

    Fastq.Support.log()

    Fastq.Support.error_if_directory(pa)

    ht = joinpath(pa, "fastp.html")

    js = joinpath(pa, "fastp.json")

    ou1 = joinpath(pa, FASTQ.TRIMMED_R1)

    ou2 = joinpath(pa, FASTQ.TRIMMED_R2)

    run(
        `fastp --detect_adapter_for_pe --thread $n_jo --json $js --html $ht --in1 $r1 --in2 $r2 --out1 $ou1 --out2 $ou2`,
    )

end

function psuedoalign(tr, n_jo, ou, r1, r2, fr, sd)

    FASTQ.Support.log()

    id = "$tr.kallisto_index"

    if !ispath(id)

        println("\nCreating kallisto index")

        run(`kallisto index --index $id $tr`)

        println("\nMade kallisto index at $id\n")

    end

    fu = ["kallisto", "quant"]

    ru = ["--threads", "$n_jo", "--index", "$id", "--output-dir", "$ou"]

    if r2 !== nothing

        println("Running paired end psuedoalignment")

        run(`$fu $ru $r1 $r2`)

    else

        println("Running single end psuedoalignment")

        run(`$fu --single --fragment-length $fr --sd $sd $ru $r1`)

    end

end

function align_cdna(al, sa, r1, r2, ge, n_jo)

    FASTQ.Support.log()

    FASTQ.Support.error_if_directory(al)

    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        println("\nMaking STAR indices, this may take a while\n")

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, `$ged`))

        end

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ged`,
        )

    end

    pr = joinpath(al, "$(sa).")

    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

end

function align_cdna_samples(ou, cd, re, n_jo; al = "transcriptome", fr = 51, sd = 0.05)

    FASTQ.Support.log()

    fq_ = FASTQ.Raw.find(cd)

    na_ = ["R1", "read1", "_1.fq"]

    naf = ""

    for fq1 in fq_

        for na in na_
            if occursin(na, fq1)

                naf = na

                nar = replace(naf, "1" => "2")

                fq2 = replace(fq1, naf => nar)

                if !isfile(fq2)

                    fq2 = nothing

                end

                sa = last(splitdir(splitext(split(fq1, naf)[1])[1]))

                println("Working on sample: $sa\n")

                pas = joinpath(ou, sa)

                if al == "transcriptome"

                    FASTQ.Raw.psuedoalign(re, n_jo, pas, fq1, fq2, fr, sd)

                elseif al == "genome"

                    FASTQ.Raw.align_cdna(pas, sa, fq1, fq2, re, n_jo)

                end

            end

        end

    end

end


function align_dna(al, sa, ba, r1, r2, ge, n_jo, me)

    FASTQ.Support.log()

    FASTQ.Support.error_if_directory(al)

    gei = "$ge.mmi"

    if !ispath(gei)

        run(`minimap2 -t $n_jo -d $gei $ge`)

    end

    tm = joinpath(al, "samtools_sort")

    du = joinpath(al, "$sa.unmarked_duplicates.bam")

    run(
        pipeline(
            `minimap2 -ax sr -t $n_jo -K $(me)G -R "@RG\tID:$sa\tSM:$sa" $gei $r1 $r2`,
            `samtools fixmate --threads $n_jo -u -m - -`,
            `samtools sort --threads $n_jo -T $tm -o $du`,
        ),
    )

    run(`samtools markdup --threads $n_jo --reference $ge --output-fmt BAM $du $ba`)

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

    rm(du)

end

function align_cdna(al, sa, r1, r2, ge, n_jo)

    FASTQ.Support.log()

    FASTQ.Support.error_if_directory(al)

    # Change to star index made with gtf gene annotation file
    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        println("\nMaking STAR indices, this may take a while\n")

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, `$ged`))

        end

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ged`,
        )

    end

    pr = joinpath(al, "$(sa).")

    # Update star run
    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

end

end
