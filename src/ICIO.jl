"""
    ICIO

Fast value-added and Global Value Chain (GVC) decompositions of Inter-Country Input-Output
(ICIO) tables, following the Borin & Mancini (2019) framework implemented by the Stata `icio`
command (Belotti, Borin & Mancini 2021).

Workflow:

```julia
using ICIO
m = read_icio_csv("EM_2015.csv", "EM_countrylist.csv")   # or load_icio(VA, FD, T; regions, sectors)
decompose(m; level = :country, perspective = :world, approach = :sink)  # 9-term corrected KWW
decompose(m; level = :sector)                                            # 13-term exporter/source
decompose(m; level = :bilateral)                                         # 13-term, per importer
```
"""
module ICIO

using LinearAlgebra
using DataFrames
import CSV
import Tables

export ICIOModel, load_icio, read_icio_csv,
       decompose, decompose_country, decompose_sector, decompose_bilateral

include("model.jl")
include("load.jl")
include("decompose.jl")
include("output.jl")

end # module ICIO
