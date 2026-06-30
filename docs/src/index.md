```@meta
CurrentModule = ICIO
```

# ICIO.jl

Fast value-added and Global Value Chain (GVC) decompositions of Inter-Country Input-Output
(ICIO) tables, following the **Borin & Mancini (2019)** framework implemented by the Stata
[`icio`](https://www.tradeconomics.com/icio/) command (Belotti, Borin & Mancini, *Stata
Journal* 2021).

The package precomputes the expensive shared matrices (notably the global Leontief inverse)
**once** per table and then computes country-, sector-, and bilateral-sector-level
decompositions vectorised across all exporters, sectors and destination pairs. A full
bilateral-sector decomposition of a 245-country × 18-sector table (≈1.08 million rows)
completes in well under a second after the one-off setup.

## Installation

```julia
using Pkg
Pkg.add("ICIO")          # once registered; until then: Pkg.develop(url = "https://github.com/SebKrantz/ICIO.jl")
```

## Quick start

```julia
using ICIO

# (a) from the icio CSV format: a headerless [T | FD] matrix + a country-list file
m = read_icio_csv("EM_2015.csv", "EM_countrylist.csv"; sectors = ["AFF", "MIN", ...])

# (b) or directly from matrices (e.g. the VA / FD / T objects of an MRIO)
m = load_icio(VA, FD, T; regions = iso3, sectors = sector_codes)

# Country level — exporter perspective, source approach, 13 terms (default)
decompose(m)

# Country level — corrected KWW / Borin-Mancini (world perspective, sink approach), 9 terms
decompose(m; perspective = :world, approach = :sink)

# Sector / bilateral level, exporter perspective; choose source or sink allocation
decompose(m; level = :sector)                          # source, 13 terms
decompose(m; level = :bilateral, approach = :sink)     # sink, 10 terms (adds VAXIM)

# Self (sectoral-bilateral) perimeter, and the importer-perspective import decomposition
decompose(m; level = :bilateral, perspective = :self)  # sectbil, 7 terms
decompose(m; flow = :imports)                          # imports by country
```

Pass a `Dict` of `year => model` to [`decompose`](@ref) to process several tables at once and
stack the results with a `:year` column.

## The decompositions

The full set of supported `flow` / `level` / `perspective` / `approach` combinations (the
columns of `icio`'s perspectives and approaches):

| `flow`      | `level`      | `perspective` / `approach`        | rows                                | terms |
|-------------|--------------|-----------------------------------|-------------------------------------|-------|
| `:exports`  | `:country`   | `:exporter` / `:source` (default) | one per exporter                    | 13    |
| `:exports`  | `:country`   | `:world` / `:source`              | one per exporter                    | 9     |
| `:exports`  | `:country`   | `:world` / `:sink`                | one per exporter                    | 9     |
| `:exports`  | `:sector`    | `:exporter` / `:source`           | one per exporter-sector             | 13    |
| `:exports`  | `:sector`    | `:exporter` / `:sink`             | one per exporter-sector             | 9     |
| `:exports`  | `:sector`    | `:self`                           | one per exporter-sector             | 7     |
| `:exports`  | `:bilateral` | `:exporter` / `:source`           | one per exporter-sector × importer  | 13    |
| `:exports`  | `:bilateral` | `:exporter` / `:sink`             | one per exporter-sector × importer  | 10    |
| `:exports`  | `:bilateral` | `:self`                           | one per exporter-sector × importer  | 7     |
| `:imports`  | `:country`   | `:importer`                       | one per importer                    | 3     |
| `:imports`  | `:bilateral` | `:importer`                       | one per (importer, VA origin)       | 2     |

The `:source` approach records value added the first time it leaves the exporter's border
(suited to production-linkage / GVC analysis); `:sink` records it the last time (suited to
final-demand analysis); the two coincide at the whole-country exporter perimeter. The `:self`
perimeter draws the boundary at the export flow itself, giving the broader Johnson (2018) /
Los et al. (2016) value-added (`DVA★ ⊇ DVAsource, DVAsink`). `:world` is country-level only.

Output is a tidy `DataFrame` of absolute values. The export term columns are
`gexp dc dva vax ref ddc fc fva fdc` (9 terms) plus `davax gvc gvcb gvcf` (13 terms) or `vaxim`
(bilateral sink), satisfying `gexp = dc + fc`, `dc = dva + ddc`, `fc = fva + fdc`,
`dva = vax + ref`, `gvc = gvcb + gvcf = gexp − davax`, and `gvcb = fc + ddc`. The import columns
are `gimp va dc` with `gimp = va + dc`.

| term | meaning |
|------|---------|
| `gexp` | gross exports |
| `dc` / `fc` | domestic / foreign content |
| `dva` / `fva` | domestic / foreign value added |
| `ddc` / `fdc` | domestic / foreign double counting |
| `vax` | domestic VA absorbed abroad (Johnson-Noguera) |
| `ref` | reflection (domestic VA returning home) |
| `davax` | domestic VA directly absorbed by the importer (source approach) |
| `vaxim` | domestic VA absorbed by the importer, incl. re-processing (sink approach; `davax ⊆ vaxim ⊆ vax`) |
| `gvc` | GVC-related trade (crosses > 1 border) |
| `gvcb` / `gvcf` | backward / forward GVC participation |
| `gimp` | gross imports (`= va + dc`) |
| `va` / `dc` | value added / double counting in imports (by origin at the bilateral level) |

## Validation

Every variant has been diffed **directly against Stata `icio`** on the EMERGING 245×18 tables and
agrees to ≈1e-6 relative (Stata's CSV output precision): world/sink and world/source (country),
exporter/source and exporter/sink (sector and bilateral, including `vaxim`), and the
importer-perspective imports. The decompositions are exactly additive (bilateral → sector →
country) and satisfy the Borin-Mancini cross-engine identities to machine precision (summed over
importers the sink DVA/FVA/VAX/REF equal the source country totals; world/source and world/sink
FVA share the same world total; imports `va + dc` = gross imports). `misc/ICIO_decomp_variants.do`
regenerates the Stata references and `misc/compare_variants_stata.jl` runs the read-only diff. The
R counterpart is the `bm()` function in the [`decompr`](https://github.com/bquast/decompr) package,
which agrees with ICIO.jl to ~1e-13.

## References

* Borin, A. & Mancini, M. (2019). *Measuring What Matters in Global Value Chains and
  Value-Added Trade.* World Bank Policy Research WP 8804.
* Belotti, F., Borin, A. & Mancini, M. (2021). *icio: Economic analysis with intercountry
  input–output tables.* The Stata Journal 21(3).
* Koopman, R., Wang, Z. & Wei, S.-J. (2014). *Tracing value-added and double counting in gross
  exports.* American Economic Review 104(2).
