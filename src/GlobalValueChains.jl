"""
    GlobalValueChains

Fast value-added and Global Value Chain (GVC) decompositions following the Borin & Mancini (2019) framework implemented by the Stata `icio`
command (Belotti, Borin & Mancini 2021).

Supports the full set of `icio` perspectives/approaches: exporter (source & sink), world
(source & sink), self (sectoral / sectoral-bilateral) export perimeters at the country, sector
and bilateral levels, plus the importer-perspective decomposition of gross imports.

Workflow:

```julia
using GlobalValueChains
m = read_icio_csv("EM_2015.csv", "EM_countrylist.csv")   # or load_icio(VA, FD, T; regions, sectors)
decompose(m)                                                       # 13-term exporter/source, by country
decompose(m; perspective = :world, approach = :sink)              # 9-term corrected KWW
decompose(m; level = :sector)                                      # 13-term exporter/source by country-sector
decompose(m; level = :bilateral, approach = :sink)                # 10-term exporter/sink (adds VAXIM), per importer
decompose(m; level = :bilateral, perspective = :self)             # sectoral-bilateral (sectbil) perimeter
decompose(m; flow = :imports)                                     # importer-perspective imports, by country
```

See [`decompose`](@ref) for the complete table of supported `flow`/`level`/`perspective`/`approach`
combinations.
"""
module GlobalValueChains

using LinearAlgebra
using DataFrames
import CSV
import Tables

export ICIOModel, load_icio, read_icio_csv,
       decompose, decompose_country, decompose_sector, decompose_bilateral, decompose_imports

include("model.jl")
include("load.jl")
include("decompose.jl")
include("decompose_world.jl")
include("decompose_sink.jl")
include("decompose_self.jl")
include("imports.jl")
include("output.jl")

end # module GlobalValueChains
