# Compare GlobalValueChains.jl bilateral-sector output against the Stata sample produced by
# misc/ICIO_decomp_bil_sample.do  ->  ${data}_GVC_BIL_SEC_SAMPLE.csv
#
# Data layout (same as misc/ICIO_decomp.jl):
#   cd to the EMERGING project root (parent of ICIO_CSV/), or set ICIO_DATA to the folder
#   that contains ${data}_2015.csv, ${data}_countrylist.csv, etc.
#
# From the EMERGING project root (default relative path):
#   julia --project=/path/to/GlobalValueChains.jl misc/ICIO_decomp_bil_sample.jl
# From anywhere (override data folder):
#   export ICIO_DATA=/path/to/EMERGING_Broad_Sectors
#   julia --project=/path/to/GlobalValueChains.jl misc/ICIO_decomp_bil_sample.jl
#
# Run misc/ICIO_decomp_bil_sample.do first to generate the Stata reference CSV.

const PKG_ROOT = dirname(@__DIR__)

import Pkg
Pkg.activate(PKG_ROOT)
using ICIO, CSV, DataFrames

# ---- paths & metadata (keep in sync with ICIO_decomp_bil_sample.do) ----
const DATA_PREFIX = "EM"
const CSV_PATH = get(ENV, "ICIO_DATA", joinpath("ICIO_CSV", "EMERGING_Broad_Sectors"))
const SAMPLE_YEAR = 2015
const SECTORS = ["AFF", "FBE", "PCM", "PSM", "TEX", "WAP", "MPR", "ELM", "TEQ", "MAN",
                 "EGW", "MIN", "SMH", "TRA", "PTE", "CON", "FIB", "PAO"]

function resolve_csv_path()
    path = abspath(CSV_PATH)
    isdir(path) || error(
        "Data directory not found: $path\n" *
        "cd to the EMERGING project root or set ICIO_DATA to the folder with " *
        "$(DATA_PREFIX)_*.csv files.")
    return path
end

base = resolve_csv_path()
clist = joinpath(base, "$(DATA_PREFIX)_countrylist.csv")
refpath = joinpath(base, "$(DATA_PREFIX)_GVC_BIL_SEC_SAMPLE.csv")

if !isfile(refpath)
    println("Reference not found: ", refpath)
    println("Run the Stata script first:  do path/to/GlobalValueChains.jl/misc/ICIO_decomp_bil_sample.do")
    exit(0)
end

ref = CSV.read(refpath, DataFrame)              # year, from_region, from_sector(idx), to_region, 13 terms
exporters = unique(ref.from_region)
importers = unique(ref.to_region)
println("Sample: exporters=", exporters, "  importers=", importers, "  rows=", nrow(ref))

m = read_icio_csv(joinpath(base, "$(DATA_PREFIX)_$(SAMPLE_YEAR).csv"), clist; sectors = SECTORS)
jl = decompose(m; level = :bilateral)
jl = jl[in.(jl.from_region, Ref(Set(exporters))) .& in.(jl.to_region, Ref(Set(importers))), :]
jl.sidx = [findfirst(==(s), SECTORS) for s in jl.from_sector]

jl_k = sort(jl, [:from_region, :sidx, :to_region])
ref_k = sort(ref, [:from_region, :from_sector, :to_region])
@assert nrow(jl_k) == nrow(ref_k) "row count mismatch: julia=$(nrow(jl_k)) stata=$(nrow(ref_k))"
@assert jl_k.from_region == ref_k.from_region && jl_k.sidx == ref_k.from_sector &&
        jl_k.to_region == ref_k.to_region "key mismatch"

cols = [:gexp, :dc, :dva, :vax, :davax, :ref, :ddc, :fc, :fva, :fdc, :gvc, :gvcb, :gvcf]
println("\nBilateral exporter/source: Julia vs Stata $(DATA_PREFIX)_GVC_BIL_SEC_SAMPLE.csv (", nrow(jl_k), " rows)")
println(rpad("col", 7), rpad("maxabs", 14), "maxrel(|val|>1)")
worst = 0.0
for c in cols
    a = Float64.(jl_k[!, c])
    b = Float64.(coalesce.(ref_k[!, c], 0.0))
    maxabs = maximum(abs.(a .- b))
    big = abs.(b) .> 1.0
    maxrel = any(big) ? maximum(abs.(a[big] .- b[big]) ./ abs.(b[big])) : 0.0
    global worst = max(worst, maxrel)
    println(rpad(string(c), 7), rpad(round(maxabs, sigdigits = 4), 14), round(maxrel, sigdigits = 3))
end
println("\nworst relative error (where |value|>1): ", round(worst, sigdigits = 3),
        worst < 1e-5 ? "   ✓ PASS" : "   ⚠ CHECK")
