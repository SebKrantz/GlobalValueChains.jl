# Validation harness for the extended Borin-Mancini decomposition variants added in GlobalValueChains.jl
# v0.2.0 (world source/sink, exporter sink, self perimeter, imports).
#
# It (1) runs every new variant on a real EMERGING table, (2) asserts the BM2019 accounting
# identities and cross-engine "anchor" equalities that pin the new engines to the already
# Stata-validated source/exporter and world/sink engines, and (3) writes CSVs (full country/
# sector variants + a bilateral sample) that can be diffed against the Stata `icio` references
# produced by misc/ICIO_decomp_variants.do.
#
# Run from the EMERGING project root (parent of ICIO_CSV/), or set ICIO_DATA:
#   julia --project=/path/to/GlobalValueChains.jl misc/ICIO_decomp_variants.jl

const PKG_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(PKG_ROOT)
using GlobalValueChains, CSV, DataFrames
using Printf

const DATA_PREFIX = "EM"
const CSV_PATH = get(ENV, "ICIO_DATA", joinpath("ICIO_CSV", "EMERGING_Broad_Sectors"))
const YEAR = 2015
const SECTORS = ["AFF", "FBE", "PCM", "PSM", "TEX", "WAP", "MPR", "ELM", "TEQ", "MAN",
                 "EGW", "MIN", "SMH", "TRA", "PTE", "CON", "FIB", "PAO"]
# bilateral sample (keep CSVs small; mirror misc/ICIO_decomp_variants.do)
const EXPORTERS = ["CHN", "DEU", "USA", "ZAF", "IND"]
const IMPORTERS = ["USA", "CHN", "DEU", "NGA"]

base = abspath(CSV_PATH)
isdir(base) || error("Data directory not found: $base (set ICIO_DATA).")
clist = joinpath(base, "$(DATA_PREFIX)_countrylist.csv")

println("Loading $(DATA_PREFIX)_$(YEAR) …"); flush(stdout)
m = read_icio_csv(joinpath(base, "$(DATA_PREFIX)_$(YEAR).csv"), clist; sectors = SECTORS)
println("  $(m.G) countries × $(m.N) sectors")

npass = Ref(0); nfail = Ref(0)
function check(name, cond)
    ok = cond === true || (cond isa Real && cond < 1e-6)
    val = cond isa Real ? @sprintf(" (%.2e)", cond) : ""
    println(rpad(ok ? "  PASS" : "  FAIL", 7), name, val)
    ok ? (npass[] += 1) : (nfail[] += 1)
    return ok
end
maxdiff(a, b) = maximum(abs.(a .- b))
byregion(df, cols) = sort(combine(groupby(df, :from_region), (cols .=> sum .=> cols)...), :from_region)

# ----- compute all variants once -----
println("Decomposing all variants …"); flush(stdout)
cty_src  = sort(decompose(m; level = :country), :country)
w_sink   = decompose(m; perspective = :world, approach = :sink)
w_src    = decompose(m; perspective = :world, approach = :source)
sec_src  = decompose(m; level = :sector, approach = :source)
sec_sink = decompose(m; level = :sector, approach = :sink)
sec_self = decompose(m; level = :sector, perspective = :self)
bil_src  = decompose(m; level = :bilateral, approach = :source)
bil_sink = decompose(m; level = :bilateral, approach = :sink)
bil_self = decompose(m; level = :bilateral, perspective = :self)
imp_cty  = decompose(m; flow = :imports)
imp_bil  = decompose(m; flow = :imports, level = :bilateral)

println("\n=== Accounting identities (real data) ===")
for (nm, df) in (("world/sink", w_sink), ("world/source", w_src),
                 ("sector/sink", sec_sink), ("sector/self", sec_self),
                 ("bilateral/sink", bil_sink), ("bilateral/self", bil_self))
    check("$nm  gexp = dc + fc", maxdiff(df.gexp, df.dc .+ df.fc))
    check("$nm  dc = dva + ddc", maxdiff(df.dc, df.dva .+ df.ddc))
    check("$nm  fc = fva + fdc", maxdiff(df.fc, df.fva .+ df.fdc))
end
check("imports  gimp = va + dc", maxdiff(imp_cty.gimp, imp_cty.va .+ imp_cty.dc))

println("\n=== Cross-engine anchors (BM2019 §5.1 / §3.2) ===")
check("world source & sink FVA equal in world total", abs(sum(w_src.fva) - sum(w_sink.fva)))
check("world FVA differ by country (not trivial)", maxdiff(w_src.fva, w_sink.fva) > 1e-3)
# sink aggregated over importers (and sectors) == source country totals
aggsec = byregion(sec_sink, [:dva, :fva, :vax, :ref])
for c in (:dva, :fva, :vax, :ref)
    check("Σ sector/sink $c == country/source $c", maxdiff(aggsec[!, c], cty_src[!, c]))
end
aggbil = byregion(bil_sink, [:dva, :fva])
check("Σ bilateral/sink dva == country/source dva", maxdiff(aggbil.dva, cty_src.dva))
check("Σ bilateral/sink fva == country/source fva", maxdiff(aggbil.fva, cty_src.fva))
# DC/FC perimeter-invariant; VAXIM nests DAVAX
check("sector/sink dc == sector/source dc", maxdiff(sec_sink.dc, sec_src.dc))
check("bilateral/sink fc == bilateral/source fc", maxdiff(bil_sink.fc, bil_src.fc))
check("DAVAX ⊆ VAXIM ⊆ VAX (bilateral/sink)",
      minimum(bil_sink.vaxim .- bil_src.davax) > -1e-6 && minimum(bil_sink.vax .- bil_sink.vaxim) > -1e-6)
# self perimeter dominance (eq. 46)
check("self dva★ ≥ source dva (bilateral)", minimum(bil_self.dva .- bil_src.dva) > -1e-6)
check("self dva★ ≥ sink dva (bilateral)", minimum(bil_self.dva .- bil_sink.dva) > -1e-6)
check("self dc == source dc (perimeter-invariant)", maxdiff(bil_self.dc, bil_src.dc))
# imports world consistency + origin additivity
check("Σ imports gimp == Σ exports", abs(sum(imp_cty.gimp) - sum(cty_src.gexp)))
aggimp = sort(combine(groupby(imp_bil, :importer), [:va, :dc] .=> sum .=> [:va, :dc]), :importer)
check("Σ origin va == importer va", maxdiff(aggimp.va, sort(imp_cty, :importer).va))

# ----- write CSVs for Stata comparison -----
println("\nWriting comparison CSVs to $base …")
keepbil(df) = filter(r -> r.from_region in EXPORTERS && r.to_region in IMPORTERS, df)
out(name, df) = (CSV.write(joinpath(base, "$(DATA_PREFIX)_$(name).csv"), df); println("  $(name).csv  ($(nrow(df)) rows)"))
out("GVC_KWW_WS_BM19", w_src)                 # country, world/source (9 terms)
out("GVC_SEC_SINK_BM19", sec_sink)            # sector, exporter/sink (9 terms)
out("GVC_SEC_SELF_BM19", sec_self)            # sector, self/sectexp (7 terms)
out("GVC_BIL_SINK_SAMPLE", keepbil(bil_sink)) # bilateral sample, exporter/sink (10 terms)
out("GVC_BIL_SELF_SAMPLE", keepbil(bil_self)) # bilateral sample, self/sectbil (7 terms)
out("GVC_IMP_BM19", imp_cty)                  # imports by country

# ----- optional: diff against Stata references if present -----
function diff_stata(name, jldf, idcols, termcols)
    ref = joinpath(base, "$(DATA_PREFIX)_$(name)_STATA.csv")
    isfile(ref) || return
    st = CSV.read(ref, DataFrame)
    j = sort(jldf, idcols); s = sort(st, idcols)
    md = maximum(maximum(abs.(Float64.(j[!, c]) .- Float64.(s[!, c]))) for c in termcols)
    check("Stata diff $name", md)
end
println("\n=== Stata reference diffs (if *_STATA.csv present) ===")
diff_stata("GVC_KWW_WS_BM19", w_src, [:country], [:gexp, :dva, :vax, :fva, :fdc])
diff_stata("GVC_SEC_SINK_BM19", sec_sink, [:from_region, :from_sector], [:gexp, :dva, :vax, :fva])
diff_stata("GVC_IMP_BM19", imp_cty, [:importer], [:gimp, :va, :dc])

println("\n", nfail[] == 0 ? "ALL $(npass[]) CHECKS PASSED" : "$(nfail[]) CHECK(S) FAILED ($(npass[]) passed)")
nfail[] == 0 || exit(1)
