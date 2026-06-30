# GlobalValueChains.jl

[![CI](https://github.com/SebKrantz/GlobalValueChains.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SebKrantz/GlobalValueChains.jl/actions/workflows/CI.yml)
[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://SebKrantz.github.io/GlobalValueChains.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Fast value-added and Global Value Chain (GVC) decompositions following the **Borin & Mancini (2019)** framework implemented by the Stata
[`icio`](https://www.tradeconomics.com/icio/) command (Belotti, Borin & Mancini, *Stata
Journal* 2021).

The package precomputes the expensive shared matrices (notably the global Leontief inverse)
**once** per table and then computes country-, sector-, and bilateral-sector-level
decompositions vectorised across all exporters, sectors and destination pairs. Where the Stata
workflow re-derives the Leontief inverse on every `icio` call — tens of thousands of calls for
a full bilateral-sector run — GlobalValueChains.jl does the full bilateral-sector decomposition of a
245-country × 18-sector table (≈1.08 million rows) in **well under a second** after the one-off
setup.

## Installation

```julia
using Pkg
Pkg.develop(path = "/path/to/GlobalValueChains.jl")   # not yet registered
```

## Quick start

```julia
using GlobalValueChains

# (a) from the icio CSV format: a headerless [T | FD] matrix + a country-list file
m = read_icio_csv("EM_2015.csv", "EM_countrylist.csv"; sectors = ["AFF","MIN", ...])

# (b) or directly from matrices (e.g. the VA / FD / T objects of an MRIO)
#     VA: length GN vector (or `nothing` to use the icio residual X .- colSums(T))
#     FD: GN×G final demand   T: GN×GN intermediate transactions
m = load_icio(VA, FD, T; regions = iso3, sectors = sector_codes)

# Country / sector / bilateral level — exporter perspective, source approach, 13 terms (default)
decompose(m)                                            # = decompose(m; level = :country)
decompose(m; level = :sector)
decompose(m; level = :bilateral)

# Country level — corrected KWW / Borin-Mancini (world perspective), sink or source approach
decompose(m; perspective = :world, approach = :sink)    # 9 terms
decompose(m; perspective = :world, approach = :source)

# Sink allocation (adds VAXIM at the bilateral level), self perimeter, and imports
decompose(m; level = :bilateral, approach = :sink)      # 10 terms
decompose(m; level = :bilateral, perspective = :self)   # sectoral-bilateral perimeter, 9 terms
decompose(m; flow = :imports)                           # importer-perspective imports
```

Convenience wrappers `decompose_country(m)`, `decompose_sector(m)`, `decompose_bilateral(m)`,
and `decompose_imports(m)` are also exported.

### Multiple years

Pass a `Dict` (label ⇒ model) to run a decomposition for several tables and stack the results
with a `:year` column — the Julia equivalent of the `foreach y in $years` loop in a Stata `.do`:

```julia
clist = "EM_countrylist.csv"
years = Dict(y => read_icio_csv("EM_$(y).csv", clist) for y in (2015, 2018, 2021, 2023))
decompose(years; level = :bilateral)
```

## The decompositions

GlobalValueChains.jl covers the full set of `icio` perspectives and approaches via the `flow`, `level`,
`perspective` and `approach` keywords:

| `flow`     | `level`     | `perspective` / `approach`        | rows                                  | terms |
|------------|-------------|-----------------------------------|---------------------------------------|-------|
| `:exports` | `:country`  | `:exporter` / `:source` (default) | one per exporter                      | 13    |
| `:exports` | `:country`  | `:world` / `:source` \| `:sink`   | one per exporter                      | 9     |
| `:exports` | `:sector`   | `:exporter` / `:source`           | one per exporter-sector               | 13    |
| `:exports` | `:sector`   | `:exporter` / `:sink`             | one per exporter-sector               | 9     |
| `:exports` | `:sector`   | `:self`                           | one per exporter-sector               | 9     |
| `:exports` | `:bilateral`| `:exporter` / `:source`           | one per exporter-sector × importer    | 13    |
| `:exports` | `:bilateral`| `:exporter` / `:sink`             | one per exporter-sector × importer    | 10    |
| `:exports` | `:bilateral`| `:self`                           | one per exporter-sector × importer    | 9     |
| `:imports` | `:country`  | `:importer`                       | one per importer                      | 3     |
| `:imports` | `:bilateral`| `:importer`                       | one per (importer, VA origin)         | 2     |

`:source` records value added the first time it leaves the exporter's border (production-linkage
view); `:sink` the last time (final-demand view; adds `vaxim` at the bilateral level); the two
coincide at the country level. `:self` draws the perimeter at the flow itself (broader Johnson
2018 / Los et al. 2016 value added). `:world` is country-level only.

Output is a tidy `DataFrame` of absolute values (same units as the table); compute shares
yourself. Identifier columns are `country` (country level) or `from_region`, `from_sector`
(and `to_region` for bilateral); imports use `importer` (and `origin`). The export term columns:

* **9 terms**: `gexp dc dva vax ref ddc fc fva fdc`
* **13 terms** (exporter/source): the above plus `davax gvc gvcb gvcf`
* **10 terms** (bilateral/sink): the 9 plus `vaxim`
* **9 terms** (self): `gexp dc dva vax ref ddc fc fva fdc`
* **imports**: `gimp va dc` (country) or `va dc` (by origin)

with the accounting identities `gexp = dc + fc`, `dc = dva + ddc`, `fc = fva + fdc`,
`dva = vax + ref`, `gvc = gvcb + gvcf = gexp − davax`, `gvcb = fc + ddc`, and `gimp = va + dc`.

| term | meaning |
|------|---------|
| `gexp` | gross exports |
| `dc` / `fc` | domestic / foreign content |
| `dva` / `fva` | domestic / foreign value added |
| `ddc` / `fdc` | domestic / foreign double counting |
| `vax` | domestic VA absorbed abroad (Johnson-Noguera) |
| `ref` | reflection (domestic VA returning home) |
| `davax` | domestic VA directly absorbed by the importer (source approach) |
| `vaxim` | domestic VA absorbed by the importer, incl. re-processing (sink; `davax ⊆ vaxim ⊆ vax`) |
| `gvc` | GVC-related trade (crosses > 1 border) |
| `gvcb` / `gvcf` | backward / forward GVC participation |
| `gimp` | gross imports (`= va + dc`) |
| `va` / `dc` | value added / double counting in imports (by VA origin at the bilateral level) |

## Validation

Most decompositions have been diffed **directly against Stata `icio`** on the EMERGING
245×18 tables and agree to ≈1e-6 relative — i.e. to Stata's CSV output precision (~7 significant
figures): world/sink and world/source at the country level, exporter/source and exporter/sink at
the sector and bilateral levels (including `vaxim`), the seven perimeter-invariant self
(sectexp/sectbil) terms, and the importer-perspective imports (`gimp`/`va`/`dc`). The self-perimeter
`vax`/`ref` (the abroad/home split of the broad self `DVA★`) are derived from the perimeter-invariant
reflection share and pass the internal identities; a direct Stata diff is generated by
`misc/compare_variants_stata.jl`. The decompositions are also exactly additive (bilateral → sector →
country) and satisfy the Borin-Mancini cross-engine identities to machine precision (summed over
importers the **sink** DVA/FVA/VAX/REF equal the **source** country totals; world/source and
world/sink FVA share the same world total; `davax ⊆ vaxim ⊆ vax`; imports `va + dc` = gross
imports). Run `Pkg.test("GlobalValueChains")` for the identity/anchor checks on a synthetic table;
`misc/ICIO_decomp_variants.{jl,do}` regenerate the Stata references, and
`misc/compare_variants_stata.jl` performs the head-to-head diff (read-only — it writes nothing).

## Method & references

Algorithmically, the source/exporter split avoids forming a separate modified Leontief inverse
per exporter: with `Mₛ = Σ_{j≠s} A_{sj} B_{js}` (an `N×N` matrix) the foreign-VA-once
coefficients are `VBfor·(I + Mₛ)⁻¹`, so the whole job is one `GN×GN` inversion plus `G` tiny
`N×N` inversions and block sums. The sink, self-perimeter and importer variants need *modified*
Leontief inverses (`B^{∤s}`, `B^{sr,n}`, `B̃^r`), but each is a low-rank change of the cached `B`,
so they reuse it via Woodbury/block updates rather than re-inverting. The world/source and
world/sink foreign VA follow Borin & Mancini (2019) eqs. (52) and (54).

* Borin, A. & Mancini, M. (2019). *Measuring What Matters in Global Value Chains and
  Value-Added Trade.* World Bank Policy Research WP 8804 (WDR 2020 background paper).
* Belotti, F., Borin, A. & Mancini, M. (2021). *icio: Economic analysis with intercountry
  input–output tables.* The Stata Journal 21(3).
* Koopman, R., Wang, Z. & Wei, S.-J. (2014). *Tracing value-added and double counting in gross
  exports.* American Economic Review 104(2).

Design influenced by the R package
[`decompr`](https://github.com/bquast/decompr) (Quast, Wang, Stolzenburg & Krantz).
