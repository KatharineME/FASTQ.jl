function align_cdna_samples(ou, cd, re, n_jo; al = "transcriptome", fr = 51, sd = 0.05)

    Fastq.support.log()

    fq_ = Fastq.fastq.find(cd)

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

                    Fastq.fastq.psuedoalign(re, n_jo, pas, fq1, fq2, fr, sd)

                elseif al == "genome"

                    Fastq.fastq.align_cdna(pas, sa, fq1, fq2, re, n_jo)

                end

            end

        end

    end


end
