# GVC decompositions of the EMERGING ICIO tables — full Julia/GlobalValueChains.jl pipeline.
#
# Julia replacement for ICIO_decomp.do. Runs all three decompositions and writes the same
# three output files. Each yearly table is read and inverted ONCE, then all three
# decompositions are computed from it (the Stata .do reloads the table for every block).
#
#   Block 1  country, world / sink        (9 terms)  -> ${data}_GVC_KWW_BM19.csv
#   Block 2  exporter-sector, exp/source  (13 terms) -> ${data}_GVC_SEC_BM19.csv
#   Block 3  bilateral-sector, exp/source (13 terms) -> ${data}_GVC_BIL_SEC_BM19.csv
#
# Data layout (same as ICIO_decomp.do):
#   cd to the EMERGING project root (parent of ICIO_CSV/), or set ICIO_DATA to the folder
#   that contains ${data}_2015.csv, ${data}_countrylist.csv, etc.
#
# From the EMERGING project root (default relative path):
#   julia --project=/path/to/GlobalValueChains.jl misc/ICIO_decomp.jl
# From anywhere (override data folder):
#   export ICIO_DATA=/path/to/EMERGING_Broad_Sectors
#   julia --project=/path/to/GlobalValueChains.jl misc/ICIO_decomp.jl

const PKG_ROOT = dirname(@__DIR__)

import Pkg
Pkg.activate(PKG_ROOT)
using GlobalValueChains, CSV, DataFrames

# ---- paths & metadata (keep in sync with ICIO_decomp.do) ----
const DATA_PREFIX = "EM"
const CSV_PATH = get(ENV, "ICIO_DATA", joinpath("ICIO_CSV", "EMERGING_Broad_Sectors"))
const YEARS = (2015, 2018, 2021, 2023)
# EMERGING broad-sector codes, in table order (sector index 1..18)
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

out_cty = joinpath(base, "$(DATA_PREFIX)_GVC_KWW_BM19.csv")
out_sec = joinpath(base, "$(DATA_PREFIX)_GVC_SEC_BM19.csv")
out_bil = joinpath(base, "$(DATA_PREFIX)_GVC_BIL_SEC_BM19.csv")
foreach(p -> isfile(p) && rm(p), (out_cty, out_sec, out_bil))

function stream!(path, df, y, first)
    insertcols!(df, 1, :year => y)
    CSV.write(path, df; append = !first, writeheader = first)
    return nrow(df)
end

for (i, y) in enumerate(YEARS)
    first = i == 1
    print("Year $y: loading…"); flush(stdout)
    t0 = time()
    m = read_icio_csv(joinpath(base, "$(DATA_PREFIX)_$(y).csv"), clist; sectors = SECTORS)

    print(" country…"); flush(stdout)
    nc = stream!(out_cty, decompose(m; level = :country, perspective = :world, approach = :sink), y, first)

    print(" sector…"); flush(stdout)
    ns = stream!(out_sec, decompose(m; level = :sector), y, first)

    print(" bilateral…"); flush(stdout)
    nb = stream!(out_bil, decompose(m; level = :bilateral), y, first)

    println("  [", round(time() - t0, digits = 1), "s]  rows: cty=$nc sec=$ns bil=$nb")
    m = nothing
    GC.gc()
end

println("Done. Wrote:")
println("  ", out_cty)
println("  ", out_sec)
println("  ", out_bil)
