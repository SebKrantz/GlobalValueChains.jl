# Tidy-DataFrame assembly and multi-year batch processing.

# country-level: one row per exporting country; term columns taken from the NamedTuple
# (9 terms for world/sink, 13 for exporter/source) in their natural order.
function _df_country(m::ICIOModel, nt::NamedTuple)
    df = DataFrame(country = m.regions)
    for (name, v) in pairs(nt)
        df[!, name] = v
    end
    return df
end

# sector-level: one row per exporting country-sector
function _df_sector(m::ICIOModel, nt::NamedTuple)
    G, N = m.G, m.N
    from_region = repeat(m.regions; inner = N)
    from_sector = repeat(m.sectors; outer = G)
    df = DataFrame(from_region = from_region, from_sector = from_sector)
    for (name, v) in pairs(nt)
        df[!, name] = v
    end
    return df
end

# bilateral-sector: one row per exporter-sector × importer (r ≠ s)
function _df_bilateral(m::ICIOModel, res::NamedTuple)
    df = DataFrame(from_region = m.regions[res.exp_g],
                   from_sector = m.sectors[res.exp_n],
                   to_region   = m.regions[res.imp_r])
    for (name, v) in pairs(res.terms)
        df[!, name] = v
    end
    return df
end

"""
    decompose(years::AbstractDict; level = :country, perspective = :exporter, approach = :source)

Batch version: `years` maps a year (or any label) to an [`ICIOModel`](@ref). Runs the
decomposition for each, prepends a `:year` column, and vertically concatenates the results
(rows ordered by sorted year). Mirrors the `foreach y in \$years` loop of `ICIO_decomp.do`.
"""
function decompose(years::AbstractDict; level::Symbol = :country,
                   perspective::Symbol = :exporter, approach::Symbol = :source)
    ks = collect(keys(years))
    try
        sort!(ks)
    catch
        # leave in iteration order if keys are not sortable
    end
    parts = DataFrame[]
    for y in ks
        df = decompose(years[y]; level = level, perspective = perspective, approach = approach)
        insertcols!(df, 1, :year => fill(y, nrow(df)))
        push!(parts, df)
    end
    return reduce(vcat, parts)
end
