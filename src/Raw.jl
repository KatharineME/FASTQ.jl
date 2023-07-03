module Raw

using FASTQ

function find(di)

    FASTQ.Support.log_sub_level_function()

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


    @info "File types found in $di"

    for na in na_n

        @info na

    end

    @info "Size of gzipped files"

    for fi in re_

        @info "File $fi is:" Base.format_bytes(stat(fi).size)

    end

    re_

end

function check_read(pa, re_, n_jo)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.error_if_directory(pa)

    th = minimum((length(re_), n_jo))

    run(`fastqc --threads $th --outdir $pa $re_`)

    run(`multiqc --outdir $pa $pa`)

end

function check_read(pa, r1, r2, n_jo; sor1 = nothing, sor2 = nothing)

    if sor1 === nothing

        re_ = [r1, r2]

    else

        re_ = [r1, r2, sor1, sor2]

    end

    check_read(joinpath(pa, "check_raw"), re_, n_jo)

end

function concatenate(fq_; na = "R1")

    FASTQ.Support.log_sub_level_function()

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

    @info "Number of forward read files = $n_fo"

    @info "Number of reverse read files = $n_re"

    sa = last(splitdir(dirname(fq_[1])))

    co = joinpath(dirname(dirname(fq_[1])), "$(sa)_concatenated")

    if n_fo <= 1 && n_re <= 1

        @warn "Nothing to concatenate, number of forward reads and reverse reads are both <= 1"

    else

        FASTQ.Support.error_if_directory(co)

        @info "Concatenating"

        gr_su = Dict(fo_ => "_R1.fastq.gz", re_ => "_R2.fastq.gz")

        for gr in keys(gr_su)

            run(pipeline(`cat $gr`, stdout = joinpath(co, "$(sa)$(gr_su[gr])")))

        end

        @info "Concatenated files saved at: $co"

    end

end

function trim(pa, n_jo, r1, r2)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.error_if_directory(pa)

    ht = joinpath(pa, "fastp.html")

    js = joinpath(pa, "fastp.json")

    ou1 = joinpath(pa, FASTQ.TR1)

    ou2 = joinpath(pa, FASTQ.TR2)

    run(
        `fastp --detect_adapter_for_pe --thread $n_jo --json $js --html $ht --in1 $r1 --in2 $r2 --out1 $ou1 --out2 $ou2`,
    )

end

function psuedoalign(ou, tr, n_jo, r1, r2, fr, sd)

    FASTQ.Support.log_sub_level_function()

    id = "$tr.kallisto_index"

    if !ispath(id)

        @info "Creating kallisto index at $id"

        run(`kallisto index --index $id $tr`)

    end

    fu_ = ["kallisto", "quant"]

    ru_ = ["--threads", "$n_jo", "--index", "$id", "--output-dir", "$ou"]

    if r2 !== nothing

        run(`$fu_ $ru_ $r1 $r2`)

    else

        run(`$fu_ --single --fragment-length $fr --sd $sd $ru_ $r1`)

    end

end

function align_cdna(pa, ge, n_jo, sa, r1, r2)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.error_if_directory(pa)

    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        @info "Making STAR indices"

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, `$ged`))

        end

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ged`,
        )

    end

    pr = joinpath(pa, "$(sa).")

    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

end

function align_cdna(pa, cd, re, n_jo; al = "transcriptome", fr = 51, sd = 0.05)

    FASTQ.Support.log_sub_level_function()

    fq_ = find(cd)

    na_ = ["R1", "read1", "_1.fq"]

    for fq1 in fq_

        for na in na_

            if occursin(na, fq1)

                fq2 = replace(fq1, na => replace(na, "1" => "2"))

                if !isfile(fq2)

                    fq2 = nothing

                end

                sa = last(splitdir(splitext(split(fq1, na)[1])[1]))

                @info "Working on sample: $sa"

                pas = joinpath(pa, sa)

                if al == "transcriptome"

                    psuedoalign(pas, re, n_jo, fq1, fq2, fr, sd)

                elseif al == "genome"

                    align_cdna(pas, re, n_jo, sa, fq1, fq2)

                end

            end

        end

    end

end


function align_dna(pa, sa, ba, r1, r2, ge, n_jo, me)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.error_if_directory(pa)

    gei = "$ge.mmi"

    if !ispath(gei)

        run(`minimap2 -t $n_jo -d $gei $ge`)

    end

    tm = joinpath(pa, "samtools_sort")

    du = joinpath(pa, "$sa.unmarked_duplicates.bam")

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

function align_single_cell_cdna(pa, sa, r1, r2, ge, n_jo)

    FASTQ.Support.log_sub_level_function()

    FASTQ.Support.error_if_directory(pa)

    # Change to star index made with gtf gene annotation file
    id = joinpath(dirname(ge), "star_indexes")

    if !ispath(id)

        mkdir(id)

        ged = splitext(ge)[1]

        if !isfile(ged)

            run(pipeline(`bgzip --decompress --stdout $ge`, `$ged`))

        end

        run(
            `star --runThreadN $n_jo --runMode genomeGenerate --genomeDir $id --genomeFastaFiles $ged`,
        )

    end

    pr = joinpath(pa, "$(sa).")

    # Update star run
    run(
        `star --runThreadN $n_jo --genomeDir $id --readFilesIn $r1 $r2 --readFilesCommand "gzip --decompress --stdout" --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $pr`,
    )

    ba = "$(pr)Aligned.sortedByCoord.out.bam"

    run(`samtools index -@ $n_jo $ba`)

    run(pipeline(`samtools stats --threads $n_jo $ba`, "$ba.stat"))

end

end
