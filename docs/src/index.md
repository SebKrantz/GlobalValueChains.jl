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

# Country level — corrected KWW / Borin-Mancini (world perspective, sink approach), 9 terms
decompose(m; level = :country, perspective = :world, approach = :sink)

# Country / sector / bilateral level — exporter perspective, source approach, 13 terms
decompose(m; level = :country)     # perspective = :exporter, approach = :source (defaults)
decompose(m; level = :sector)
decompose(m; level = :bilateral)
```

Pass a `Dict` of `year => model` to [`decompose`](@ref) to process several tables at once and
stack the results with a `:year` column.

## The decompositions

| `level`      | `perspective` / `approach`        | rows                                | terms |
|--------------|-----------------------------------|-------------------------------------|-------|
| `:country`   | `:world` / `:sink`                | one per exporter                    | 9     |
| `:country`   | `:exporter` / `:source` (default) | one per exporter                    | 13    |
| `:sector`    | `:exporter` / `:source`           | one per exporter-sector             | 13    |
| `:bilateral` | `:exporter` / `:source`           | one per exporter-sector × importer  | 13    |

Output is a tidy `DataFrame` of absolute values. The term columns are
`gexp dc dva vax ref ddc fc fva fdc` (9 terms) plus `davax gvc gvcb gvcf` (13 terms), satisfying
`gexp = dc + fc`, `dc = dva + ddc`, `fc = fva + fdc`, `dva = vax + ref`,
`gvc = gvcb + gvcf = gexp − davax`, and `gvcb = fc + ddc`.

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

The country-level world/sink decomposition reproduces the Stata `icio` output to ~1e-7 relative
error; the source/exporter decomposition satisfies all the accounting identities above and is
exactly additive (bilateral → sector → country). The R counterpart is the `bm()` function in the
[`decompr`](https://github.com/bquast/decompr) package, which agrees with ICIO.jl to ~1e-13.

## References

* Borin, A. & Mancini, M. (2019). *Measuring What Matters in Global Value Chains and
  Value-Added Trade.* World Bank Policy Research WP 8804.
* Belotti, F., Borin, A. & Mancini, M. (2021). *icio: Economic analysis with intercountry
  input–output tables.* The Stata Journal 21(3).
* Koopman, R., Wang, Z. & Wei, S.-J. (2014). *Tracing value-added and double counting in gross
  exports.* American Economic Review 104(2).
