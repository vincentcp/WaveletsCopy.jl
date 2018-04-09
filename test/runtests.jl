using WaveletsCopy
using Base.Test
CREATE_README = true
include("suite_dwtstep.jl")
include("suite_evaluation.jl")

if CREATE_README
    try
        println("Create README.md")
        run(`jupyter nbconvert --execute --to markdown --output README.md notebooks/README.ipynb`)
        run(`mv notebooks/README.md .`)
        try run(`rm -rf README_files/`) end
        run(`mv notebooks/README_files/ README_files/`)
    catch
        nothing
    end
end

print("\ntesting: success\n")
