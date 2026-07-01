# READ-ONLY comparison of GlobalValueChains.jl against the Stata `icio` reference files produced by
# misc/ICIO_decomp_variants.do. Computes each Julia decomposition in memory and diffs it
# against the matching EM_GVC_*_STATA.csv. It WRITES NOTHING — no output files are created or
# overwritten, so it cannot clobber the Stata references.
#
# Run from the EMERGING project root (or set ICIO_DATA):
#   julia --project=/path/to/GlobalValueChains.jl misc/compare_variants_stata.jl

const PKG_ROOT = dirname(@__DIR__)
import Pkg; Pkg.activate(PKG_ROOT)
using GlobalValueChains, CSV, DataFrames

const DATA_PREFIX = "EM"
const CSV_PATH = get(ENV, "ICIO_DATA", joinpath("ICIO_CSV", "EMERGING_Broad_Sectors"))
const YEAR = 2015
const SECTORS = ["AFF", "FBE", "PCM", "PSM", "TEX", "WAP", "MPR", "ELM", "TEQ", "MAN",
                 "EGW", "MIN", "SMH", "TRA", "PTE", "CON", "FIB", "PAO"]
const SECIDX = Dict(s => i for (i, s) in enumerate(SECTORS))
const EXPORTERS = ["CHN", "DEU", "USA", "ZAF", "IND"]
const IMPORTERS = ["USA", "CHN", "DEU", "NGA"]

base = abspath(CSV_PATH)
m = read_icio_csv(joinpath(base, "$(DATA_PREFIX)_$(YEAR).csv"),
                  joinpath(base, "$(DATA_PREFIX)_countrylist.csv"); sectors = SECTORS)
println("Loaded $(DATA_PREFIX)_$(YEAR): $(m.G)×$(m.N).  (this script writes no files)\n")

# diff one variant; `jl` already restricted/keyed to match the Stata layout
function compare(label, file, jl, idcols, terms)
    path = joinpath(base, file)
    if !isfile(path)
        println("• $label: $file not found — run misc/ICIO_decomp_variants.do first."); return
    end
    st = CSV.read(path, DataFrame)
    if nrow(st) == 0 || all(ismissing, st[!, terms[1]])
        println("• $label: $file is empty — re-run the (fixed) .do to populate it."); return
    end
    rename!(st, Dict(t => Symbol(t, "_s") for t in terms))
    j = innerjoin(jl, st, on = idcols)
    println("• $label  ($(nrow(j))/$(nrow(st)) rows matched)")
    for t in terms
        a = Float64.(j[!, t]); b = Float64.(j[!, Symbol(t, "_s")])
        mx = maximum(abs.(a .- b)); rel = maximum(abs.(a .- b) ./ (abs.(b) .+ 1))
        println("    ", rpad(t, 6), "max|Δ|=", rpad(round(mx, sigdigits = 3), 11),
                "max relΔ=", round(rel, sigdigits = 3))
    end
end

# 1) country, world / source
compare("country world/source", "$(DATA_PREFIX)_GVC_KWW_WS_BM19_STATA.csv",
        decompose(m; perspective = :world, approach = :source),
        [:country], [:gexp, :dc, :dva, :vax, :ref, :ddc, :fc, :fva, :fdc])

# 2) sector, exporter / sink
sec = decompose(m; level = :sector, approach = :sink)
sec.from_sector = [SECIDX[s] for s in sec.from_sector]
compare("sector exporter/sink", "$(DATA_PREFIX)_GVC_SEC_SINK_BM19_STATA.csv",
        sec, [:from_region, :from_sector], [:gexp, :dc, :dva, :vax, :ref, :ddc, :fc, :fva, :fdc])

# 3) bilateral, exporter / sink (sample) — includes VAXIM
bil = filter(r -> r.from_region in EXPORTERS && r.to_region in IMPORTERS,
             decompose(m; level = :bilateral, approach = :sink))
bil.from_sector = [SECIDX[s] for s in bil.from_sector]
compare("bilateral exporter/sink", "$(DATA_PREFIX)_GVC_BIL_SINK_SAMPLE_STATA.csv",
        bil, [:from_region, :from_sector, :to_region],
        [:gexp, :dc, :dva, :vax, :vaxim, :ref, :ddc, :fc, :fva, :fdc])

# 4) imports, importer perspective
compare("imports importer", "$(DATA_PREFIX)_GVC_IMP_BM19_STATA.csv",
        decompose(m; flow = :imports), [:importer], [:gimp, :va, :dc])

# 5) sector, self (sectexp) perimeter — 9 terms
secs = decompose(m; level = :sector, perspective = :self)
secs.from_sector = [SECIDX[s] for s in secs.from_sector]
compare("sector self (sectexp)", "$(DATA_PREFIX)_GVC_SEC_SELF_SAMPLE_STATA.csv",
        secs, [:from_region, :from_sector], [:gexp, :dc, :dva, :vax, :ref, :ddc, :fc, :fva, :fdc])

# 6) bilateral, self (sectbil) perimeter (sample) — 9 terms
bils = filter(r -> r.from_region in EXPORTERS && r.to_region in IMPORTERS,
              decompose(m; level = :bilateral, perspective = :self))
bils.from_sector = [SECIDX[s] for s in bils.from_sector]
compare("bilateral self (sectbil)", "$(DATA_PREFIX)_GVC_BIL_SELF_SAMPLE_STATA.csv",
        bils, [:from_region, :from_sector, :to_region],
        [:gexp, :dc, :dva, :vax, :ref, :ddc, :fc, :fva, :fdc])
