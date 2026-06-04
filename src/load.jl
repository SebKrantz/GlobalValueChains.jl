# Reading the icio CSV table format.

"""
    read_icio_csv(table_csv, countrylist_csv; sectors = nothing, X = nothing, VA = nothing)

Load an [`ICIOModel`](@ref) from the icio CSV format used by the Stata `icio` command and
produced by `STATA_ICIO_CSVs_V2.R`:

* `table_csv` — headerless `GN × (GN + G)` matrix `[T | FD]`: the first `GN` columns are the
  intermediate transactions `T`, the last `G` columns the final demand `FD` (one per country).
* `countrylist_csv` — headerless one-column file of `G` country codes (e.g. ISO3).

The number of sectors is inferred as `N = GN / G`. `sectors` names the industries and may be a
vector of `N` codes, a path to a headerless one-column CSV of sector codes (like
`countrylist_csv`), or `nothing` (defaults to `"sector1"…"sectorN"`). Supplying the real sector
codes means every output `DataFrame` carries them in `from_sector` from the start. By default
value added is the icio column residual (`X .- colSums(T)`), reproducing `icio` exactly; pass
`VA`/`X` to override.
"""
function read_icio_csv(table_csv::AbstractString, countrylist_csv::AbstractString;
                       sectors::Union{Nothing,AbstractVector,AbstractString} = nothing,
                       X::Union{Nothing,AbstractVector} = nothing,
                       VA::Union{Nothing,AbstractVector} = nothing)
    regions = String.(Tables.getcolumn(CSV.File(countrylist_csv; header = false), 1))
    G = length(regions)

    M = Tables.matrix(CSV.File(table_csv; header = false, types = Float64))
    GN, ncol = size(M)
    GN % G == 0 || error("Table has $GN rows but $G countries; $GN is not divisible by $G.")
    N = GN ÷ G
    ncol == GN + G || error(
        "Table has $ncol columns, expected GN + G = $(GN + G) (GN intermediate + G final demand).")

    T = M[:, 1:GN]
    FD = M[:, GN+1:GN+G]

    secs = if sectors === nothing
        ["sector$n" for n in 1:N]
    elseif sectors isa AbstractString
        String.(Tables.getcolumn(CSV.File(sectors; header = false), 1))
    else
        collect(String, sectors)
    end
    length(secs) == N || throw(DimensionMismatch(
        "got $(length(secs)) sector names, expected N = $N"))

    return load_icio(VA, FD, T; regions = regions, sectors = secs, X = X)
end
