pas = "/Users/kate/Downloads/ferreira_treg/output/7166-MR-100/"

# Use the raw STARsolo output. If false, the filtered output will be used.
ra = true

or = "human"

using BioLab

using CSV

using DataFrames

function make_gene_by_sample(pas, fi, ma)

    if ra == true

        pa = "Solo.out/Gene/raw/"

    else

        pa = "Solo.out/Gene/filtered/"

    end

    # ma = CSV.read(joinpath(pas, pa, "matrix.mtx"), DataFrame; delim='\t')

    # ce = CSV.read(joinpath(pas, pa, "barcodes.tsv"), DataFrame; delim='\t')
    
    # ge = CSV.read(joinpath(pas, pa, "features.tsv"), DataFrame; delim='\t')
    
    ma = DataFrame(feature=[1, 5, 3], cell=[1, 1, 2], count=[1, 15, 7])

    ge_ = ["A", "B", "C", "D", "E"]

    ce_ = ["cell1", "cell2", "cell3"]

    df = DataFrame()

    df.gene = 


    ma
    # Build gene by cell
    
    # Plot expression per gene conveniently 
    
    # Plot mitochondrial expression conveniently
    
end
