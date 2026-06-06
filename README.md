# ICIO.jl

[![CI](https://github.com/SebKrantz/ICIO.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SebKrantz/ICIO.jl/actions/workflows/CI.yml)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://SebKrantz.github.io/ICIO.jl/dev/)

Fast value-added and Global Value Chain (GVC) decompositions of Inter-Country Input-Output
(ICIO) tables, following the **Borin & Mancini (2019)** framework implemented by the Stata
[`icio`](https://www.tradeconomics.com/icio/) command (Belotti, Borin & Mancini, *Stata
Journal* 2021).

The package precomputes the expensive shared matrices (notably the global Leontief inverse)
**once** per table and then computes country-, sector-, and bilateral-sector-level
decompositions vectorised across all exporters, sectors and destination pairs. Where the Stata
workflow re-derives the Leontief inverse on every `icio` call — tens of thousands of calls for
a full bilateral-sector run — ICIO.jl does the full bilateral-sector decomposition of a
245-country × 18-sector table (≈1.08 million rows) in **well under a second** after the one-off
setup.

## Installation

```julia
using Pkg
Pkg.develop(path = "/path/to/ICIO.jl")   # not yet registered
```

## Quick start

```julia
using ICIO

# (a) from the icio CSV format: a headerless [T | FD] matrix + a country-list file
m = read_icio_csv("EM_2015.csv", "EM_countrylist.csv"; sectors = ["AFF","MIN", ...])

# (b) or directly from matrices (e.g. the VA / FD / T objects of an MRIO)
#     VA: length GN vector (or `nothing` to use the icio residual X .- colSums(T))
#     FD: GN×G final demand   T: GN×GN intermediate transactions
m = load_icio(VA, FD, T; regions = iso3, sectors = sector_codes)

# Country level — corrected KWW / Borin-Mancini (world perspective, sink approach), 9 terms
decompose(m; level = :country, perspective = :world, approach = :sink)

# Country / sector / bilateral level — exporter perspective, source approach, 13 terms
decompose(m; level = :country)     # perspective = :exporter, approach = :source (defaults)
decompose(m; level = :sector)
decompose(m; level = :bilateral)
```

Convenience wrappers `decompose_country(m)`, `decompose_sector(m)`, `decompose_bilateral(m)`
are also exported.

### Multiple years

Pass a `Dict` (label ⇒ model) to run a decomposition for several tables and stack the results
with a `:year` column — the Julia equivalent of the `foreach y in $years` loop in a Stata `.do`:

```julia
clist = "EM_countrylist.csv"
years = Dict(y => read_icio_csv("EM_$(y).csv", clist) for y in (2015, 2018, 2021, 2023))
decompose(years; level = :bilateral)
```

## The decompositions

| `level`     | `perspective` / `approach`        | rows                                  | terms |
|-------------|-----------------------------------|---------------------------------------|-------|
| `:country`  | `:world` / `:sink`                | one per exporter                      | 9     |
| `:country`  | `:exporter` / `:source` (default) | one per exporter                      | 13    |
| `:sector`   | `:exporter` / `:source`           | one per exporter-sector               | 13    |
| `:bilateral`| `:exporter` / `:source`           | one per exporter-sector × importer    | 13    |

Output is a tidy `DataFrame` of absolute values (same units as the table); compute shares
yourself. Identifier columns are `country` (country level) or `from_region`, `from_sector`
(and `to_region` for bilateral). The term columns are:

* **9 terms** (world/sink): `gexp dc dva vax ref ddc fc fva fdc`
* **13 terms** (exporter/source): the above plus `davax gvc gvcb gvcf`

with the accounting identities `gexp = dc + fc`, `dc = dva + ddc`, `fc = fva + fdc`,
`dva = vax + ref`, `gvc = gvcb + gvcf = gexp − davax`, and `gvcb = fc + ddc`.

| term | meaning |
|------|---------|
| `gexp` | gross exports |
| `dc` / `fc` | domestic / foreign content |
| `dva` / `fva` | domestic / foreign value added |
| `ddc` / `fdc` | domestic / foreign double counting |
| `vax` | domestic VA absorbed abroad (Johnson-Noguera) |
| `ref` | reflection (domestic VA returning home) |
| `davax` | domestic VA directly absorbed by the importer |
| `gvc` | GVC-related trade (crosses > 1 border) |
| `gvcb` / `gvcf` | backward / forward GVC participation |

## Validation

The country-level world/sink decomposition reproduces the Stata `icio` output
(`EM_GVC_KWW_BM19.csv`) to ~1e-7 relative error on all nine terms. The source/exporter
decomposition satisfies all the accounting identities above and is exactly additive
(bilateral → sector → country) to floating-point precision; its perspective-invariant terms
(`gexp dc dva vax ref ddc fc`) match the Stata reference. Run `Pkg.test("ICIO")` for the
identity and additivity checks on a small synthetic table.

## Method & references

Algorithmically, the source/exporter split avoids forming a separate modified Leontief inverse
per exporter: with `Mₛ = Σ_{j≠s} A_{sj} B_{js}` (an `N×N` matrix) the foreign-VA-once
coefficients are `VBfor·(I + Mₛ)⁻¹`, so the whole job is one `GN×GN` inversion plus `G` tiny
`N×N` inversions and block sums. The world/sink foreign VA follows Borin & Mancini (2019)
eq. (54).

* Borin, A. & Mancini, M. (2019). *Measuring What Matters in Global Value Chains and
  Value-Added Trade.* World Bank Policy Research WP 8804 (WDR 2020 background paper).
* Belotti, F., Borin, A. & Mancini, M. (2021). *icio: Economic analysis with intercountry
  input–output tables.* The Stata Journal 21(3).
* Koopman, R., Wang, Z. & Wei, S.-J. (2014). *Tracing value-added and double counting in gross
  exports.* American Economic Review 104(2).

Design influenced by the R package
[`decompr`](https://github.com/bquast/decompr) (Quast, Wang, Stolzenburg & Krantz).
