using WaveletsCopy

CREATE_README = false
include("suite_dwtstep.jl")
include("suite_evaluation.jl")
include("suite_quadrature.jl")

if CREATE_README
    try
        println("Create README.md")
        run(`jupyter nbconvert --execute --to markdown --output README.md notebooks/README.ipynb`)
        run(`mv notebooks/README.md .`)
        try
            run(`rm -rf README_files/`)
        finally
            run(`mv notebooks/README_files/ README_files/`)
        end
    catch
        nothing
    end
end

print("\ntesting: success\n")
