# Value-added / GVC decompositions.
#
# Two engines:
#  * source/exporter perspective (1b country, 2 sector, 3 bilateral) — 13 terms,
#    additive across destinations and sectors. Uses cached coefficients + cheap per-pair work.
#  * world/sink perspective (1a country) — 9 terms, a vectorised port of decompr::kww
#    (terms T1..T9 of Koopman-Wang-Wei, corrected), aggregated to the country level.

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
_BFD(m::ICIOModel) = m.B * m.FD   # output driven by each absorbing country's final demand

# Wcol[r-block] = L_rr * FD[r-block, r]  (local output of r for its own final demand)
function _Wcol(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    Wcol = zeros(GN)
    for r in 1:G
        rng = blockrange(r, N)
        @view(Wcol[rng]) .= @view(m.L[rng, rng]) * @view(m.FD[rng, r])
    end
    return Wcol
end

# ---------------------------------------------------------------------------
# Source / exporter perspective, sector level (13 terms per country-sector)
# ---------------------------------------------------------------------------
function _source_sector(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, L, FD, X, E, ESR = m.A, m.L, m.FD, m.X, m.E, m.ESR
    VBdom, VBfor, VLdom, fvacoef = m.VBdom, m.VBfor, m.VLdom, m.fvacoef
    BFD = _BFD(m)
    Wcol = _Wcol(m)

    DAEsum = zeros(GN)   # Σ_{r≠s} DAE_sr   (per export country-sector)
    VAXEsum = zeros(GN)  # Σ_{r≠s} VAXE_sr
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        for r in 1:G
            r == s && continue
            rrng = blockrange(r, N)
            Asr = @view A[srng, rrng]
            dae = Asr * @view(Wcol[rrng])                       # N-vector
            xnotabs = @view(X[rrng]) .- @view(BFD[rrng, s])     # N-vector
            vaxe = Asr * xnotabs                                # N-vector
            for (k, n) in enumerate(srng)
                DAEsum[n]  += FD[n, r] + dae[k]
                VAXEsum[n] += FD[n, r] + vaxe[k]
            end
        end
    end

    gexp = copy(E)
    dc  = VBdom .* E
    fc  = VBfor .* E
    dva = VLdom .* E
    ddc = dc .- dva
    fva = fvacoef .* E
    fdc = fc .- fva
    davax = VLdom .* DAEsum
    vax = VLdom .* VAXEsum
    ref = dva .- vax
    gvc = E .- davax
    gvcb = fc .+ ddc
    gvcf = gvc .- gvcb

    return (gexp = gexp, dc = dc, dva = dva, vax = vax, davax = davax, ref = ref,
            ddc = ddc, fc = fc, fva = fva, fdc = fdc, gvc = gvc, gvcb = gvcb, gvcf = gvcf)
end

# country-level source decomposition = sector terms summed within country
function _source_country(m::ICIOModel)
    s = _source_sector(m)
    G, N = m.G, m.N
    agg(v) = (out = zeros(G); @inbounds for i in eachindex(v); out[ctry(i, N)] += v[i]; end; out)
    return map(agg, s)
end

# ---------------------------------------------------------------------------
# Source / exporter perspective, bilateral-sector level (13 terms per s,n,r with r≠s)
# Returns index vectors (exporter, sector, importer) plus the 13 term vectors.
# ---------------------------------------------------------------------------
function _source_bilateral(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, L, FD, X, ESR = m.A, m.L, m.FD, m.X, m.ESR
    VBdom, VBfor, VLdom, fvacoef = m.VBdom, m.VBfor, m.VLdom, m.fvacoef
    BFD = _BFD(m)
    Wcol = _Wcol(m)

    nrow = G * (G - 1) * N
    exp_g = Vector{Int}(undef, nrow)
    exp_n = Vector{Int}(undef, nrow)
    imp_r = Vector{Int}(undef, nrow)
    gexp = Vector{Float64}(undef, nrow); dc = similar(gexp); dva = similar(gexp)
    vax = similar(gexp); davax = similar(gexp); ref = similar(gexp); ddc = similar(gexp)
    fc = similar(gexp); fva = similar(gexp); fdc = similar(gexp)
    gvc = similar(gexp); gvcb = similar(gexp); gvcf = similar(gexp)

    row = 0
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        for r in 1:G
            r == s && continue
            rrng = blockrange(r, N)
            Asr = @view A[srng, rrng]
            dae = Asr * @view(Wcol[rrng])
            xnotabs = @view(X[rrng]) .- @view(BFD[rrng, s])
            vaxe = Asr * xnotabs
            for (k, n) in enumerate(srng)
                row += 1
                e = ESR[n, r]
                DCv = VBdom[n] * e; FCv = VBfor[n] * e
                DVAv = VLdom[n] * e; DDCv = DCv - DVAv
                FVAv = fvacoef[n] * e
                DAVAXv = VLdom[n] * (FD[n, r] + dae[k])
                VAXv = VLdom[n] * (FD[n, r] + vaxe[k])
                GVCv = e - DAVAXv; GVCBv = FCv + DDCv
                exp_g[row] = s; exp_n[row] = k; imp_r[row] = r
                gexp[row] = e
                dc[row] = DCv; fc[row] = FCv
                dva[row] = DVAv; ddc[row] = DDCv
                fva[row] = FVAv; fdc[row] = FCv - FVAv
                davax[row] = DAVAXv; vax[row] = VAXv; ref[row] = DVAv - VAXv
                gvc[row] = GVCv; gvcb[row] = GVCBv; gvcf[row] = GVCv - GVCBv
            end
        end
    end

    return (exp_g = exp_g, exp_n = exp_n, imp_r = imp_r,
            terms = (gexp = gexp, dc = dc, dva = dva, vax = vax, davax = davax, ref = ref,
                     ddc = ddc, fc = fc, fva = fva, fdc = fdc, gvc = gvc, gvcb = gvcb, gvcf = gvcf))
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
"""
    decompose(m::ICIOModel; flow = :exports, level = :country,
              perspective = :exporter, approach = :source)

Decompose gross trade into value-added / GVC components (Borin & Mancini 2019), returning a
tidy `DataFrame`. Mirrors the Stata `icio` command's perspectives and approaches.

`flow = :exports` (default) decomposes exports:

| `level`     | `perspective` | `approach`      | terms | description |
|:------------|:--------------|:----------------|:------|:------------|
| `:country`  | `:exporter`   | `:source`(=sink)| 13    | gexp dc dva vax davax ref ddc fc fva fdc gvc gvcb gvcf |
| `:country`  | `:world`      | `:source`       | 9     | world perimeter, FVA at first foreign crossing (eq. 52) |
| `:country`  | `:world`      | `:sink`         | 9     | world perimeter, corrected KWW (eq. 54) |
| `:sector`   | `:exporter`   | `:source`       | 13    | country perimeter, sectoral breakdown |
| `:sector`   | `:exporter`   | `:sink`         | 9     | gexp dc dva vax ref ddc fc fva fdc |
| `:sector`   | `:self`       | —               | 9     | sectoral (sectexp) perimeter |
| `:bilateral`| `:exporter`   | `:source`       | 13    | one row per exporter-sector × importer (r≠s) |
| `:bilateral`| `:exporter`   | `:sink`         | 10    | adds `vaxim` (DVA absorbed by direct importer, eq. 39) |
| `:bilateral`| `:self`       | —               | 9     | sectoral-bilateral (sectbil) perimeter |

The `:source` approach records value added the first time it leaves country `s`'s border (suited
to production-linkage / GVC analysis); `:sink` records it the last time (suited to final-demand
analysis). At the whole-country exporter perimeter the two coincide. The `:self` perimeter uses
the broader Johnson (2018) / Los et al. (2016) value-added notion (`DVA★ ⊇ DVAsource, DVAsink`).
`:world` is available at the country level only.

`flow = :imports` decomposes a country's gross imports from the importer perspective (eq. 51):

| `level`     | terms        | description |
|:------------|:-------------|:------------|
| `:country`  | gimp va dc   | one row per importer |
| `:bilateral`| va dc        | one row per (importer, value-added origin); sums over origin to imports |

Columns are absolute values (same currency units as the table). See [`read_icio_csv`](@ref) /
[`load_icio`](@ref) to construct `m`, and [`decompose`](@ref) over a `Dict` of years for batches.
"""
function decompose(m::ICIOModel; flow::Symbol = :exports, level::Symbol = :country,
                   perspective::Symbol = :exporter, approach::Symbol = :source)
    flow === :imports && return _decompose_imports(m, level, perspective)
    flow === :exports || error("flow must be :exports or :imports; got :$flow.")
    if level === :country
        if perspective === :world
            (approach === :sink || approach === :source) ||
                error("level=:country, perspective=:world requires approach=:sink or :source.")
            return _df_country(m, _world_country(m, approach))
        elseif perspective === :exporter
            return _df_country(m, _source_country(m))  # source ≡ sink at the country perimeter
        else
            error("flow=:exports, level=:country supports perspective=:exporter or :world; " *
                  "got perspective=:$perspective (:self/:world-bilateral are not country-level).")
        end
    elseif level === :sector
        if perspective === :exporter
            approach === :source && return _df_sector(m, _source_sector(m))
            approach === :sink   && return _df_sector(m, _sink_sector(m))
            error("approach must be :source or :sink; got :$approach.")
        elseif perspective === :self
            return _df_sector(m, _self_sector(m))
        else
            error("flow=:exports, level=:sector supports perspective=:exporter " *
                  "(approach :source/:sink) or :self; :world is country-only.")
        end
    elseif level === :bilateral
        if perspective === :exporter
            approach === :source && return _df_bilateral(m, _source_bilateral(m))
            approach === :sink   && return _df_bilateral(m, _sink_bilateral(m))
            error("approach must be :source or :sink; got :$approach.")
        elseif perspective === :self
            return _df_bilateral(m, _self_bilateral(m))
        else
            error("flow=:exports, level=:bilateral supports perspective=:exporter " *
                  "(approach :source/:sink) or :self; :world is country-only.")
        end
    else
        error("level must be :country, :sector or :bilateral; got :$level.")
    end
end

# imports router (perspective :importer; approach not applicable)
function _decompose_imports(m::ICIOModel, level::Symbol, perspective::Symbol)
    perspective === :self &&
        error("flow=:imports, sectoral-importer (sectimp) perspective is not yet implemented.")
    # :exporter is the global default; for imports the only perimeter is the importer's border
    (perspective === :importer || perspective === :exporter) ||
        error("flow=:imports supports perspective=:importer; got perspective=:$perspective.")
    if level === :country
        return _df_imports_country(m, _imports_country(m))
    elseif level === :bilateral
        return _df_imports_bilateral(m, _imports_bilateral(m))
    elseif level === :sector
        error("flow=:imports, level=:sector (sectoral imports) is not yet implemented; " *
              "use level=:country or :bilateral.")
    else
        error("flow=:imports supports level=:country or :bilateral; got :$level.")
    end
end

"Convenience wrapper for [`decompose`](@ref) at `level=:country`."
decompose_country(m::ICIOModel; perspective = :world, approach = :sink) =
    decompose(m; level = :country, perspective = perspective, approach = approach)
"Convenience wrapper for [`decompose`](@ref) at `level=:sector`."
decompose_sector(m::ICIOModel; perspective = :exporter, approach = :source) =
    decompose(m; level = :sector, perspective = perspective, approach = approach)
"Convenience wrapper for [`decompose`](@ref) at `level=:bilateral`."
decompose_bilateral(m::ICIOModel; perspective = :exporter, approach = :source) =
    decompose(m; level = :bilateral, perspective = perspective, approach = approach)
"Convenience wrapper for the import decomposition (`flow=:imports`)."
decompose_imports(m::ICIOModel; level = :country) =
    decompose(m; flow = :imports, level = level, perspective = :importer)
