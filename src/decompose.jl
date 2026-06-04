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
# World / sink foreign value added (Borin-Mancini 2019 eq. 54), country level.
# FVA_sr = VBfor·DAE_sr (foreign VA directly absorbed by importer r)
#        + V_r B_rs · (A_sr L_rr REX_r)  (importer r's own VA that r re-exports onward).
# Summed over destinations r ≠ s. FDC = FC − FVA.
# ---------------------------------------------------------------------------
function _fva_world_sink(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, B, L, Vc, FD, VBfor = m.A, m.B, m.L, m.V, m.FD, m.VBfor
    Wcol = _Wcol(m)

    # Wrex_r = L_rr * REX_r,  REX_r = Σ_{j≠r} (FD[r,j] + A_rj Wcol[j])  (r's directly-absorbed re-exports)
    Wrex = zeros(GN)
    for r in 1:G
        rrng = blockrange(r, N)
        rex = zeros(N)
        for j in 1:G
            j == r && continue
            jrng = blockrange(j, N)
            rex .+= @view(FD[rrng, j]) .+ @view(A[rrng, jrng]) * @view(Wcol[jrng])
        end
        @view(Wrex[rrng]) .= @view(L[rrng, rrng]) * rex
    end

    fva = zeros(G)
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        acc = 0.0
        for r in 1:G
            r == s && continue
            rrng = blockrange(r, N)
            Asr = @view A[srng, rrng]
            dae = @view(FD[srng, r]) .+ Asr * @view(Wcol[rrng])      # directly-absorbed exports s→r
            # term A+B: foreign VA multiplier · directly-absorbed exports
            for (k, n) in enumerate(srng)
                acc += VBfor[n] * dae[k]
            end
            # term C: importer r's own VA in s's exports that r re-exports onward
            AsrWrex = Asr * @view(Wrex[rrng])                        # N-vector (s-sectors)
            vbr = transpose(@view(B[rrng, srng])) * @view(Vc[rrng])  # N-vector (s-sectors): Σ_{m∈r} V_m B[m, s·]
            for k in 1:N
                acc += vbr[k] * AsrWrex[k]
            end
        end
        fva[s] = acc
    end
    return fva
end

# World / sink country decomposition (9 terms). Domestic side is perspective-invariant
# (taken from the source decomposition); only FVA/FDC use the world/sink split.
function _world_sink_country(m::ICIOModel)
    sc = _source_country(m)   # gexp, dc, dva, vax, davax, ref, ddc, fc, fva, fdc, gvc, gvcb, gvcf
    fva = _fva_world_sink(m)
    fdc = sc.fc .- fva
    return (gexp = sc.gexp, dc = sc.dc, dva = sc.dva, vax = sc.vax, ref = sc.ref,
            ddc = sc.ddc, fc = sc.fc, fva = fva, fdc = fdc)
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
    decompose(m::ICIOModel; level = :country, perspective = :exporter, approach = :source)

Decompose gross exports into value-added / GVC components, returning a tidy `DataFrame`.

* `level = :country`  — one row per exporting country.
  - `perspective = :world,  approach = :sink`   → 9 terms (corrected KWW / Borin-Mancini).
  - `perspective = :exporter, approach = :source` (default) → 13 terms.
* `level = :sector`   — one row per exporting country-sector (`:exporter`/`:source`, 13 terms).
* `level = :bilateral`— one row per exporter-sector × importer (`:exporter`/`:source`, 13 terms),
  excluding own-country destinations.

Columns are absolute values (same currency units as the table). See [`read_icio_csv`](@ref)
and [`load_icio`](@ref) for constructing `m`, and [`decompose`](@ref) over a `Dict` of years
for batch processing.
"""
function decompose(m::ICIOModel; level::Symbol = :country,
                   perspective::Symbol = :exporter, approach::Symbol = :source)
    if level === :country
        if perspective === :world && approach === :sink
            return _df_country(m, _world_sink_country(m))
        elseif perspective === :exporter && approach === :source
            return _df_country(m, _source_country(m))
        else
            error("level=:country supports (perspective=:world, approach=:sink) or " *
                  "(perspective=:exporter, approach=:source); got " *
                  "(perspective=:$perspective, approach=:$approach).")
        end
    elseif level === :sector
        (perspective === :exporter && approach === :source) ||
            error("level=:sector requires perspective=:exporter, approach=:source.")
        return _df_sector(m, _source_sector(m))
    elseif level === :bilateral
        (perspective === :exporter && approach === :source) ||
            error("level=:bilateral requires perspective=:exporter, approach=:source.")
        return _df_bilateral(m, _source_bilateral(m))
    else
        error("level must be :country, :sector or :bilateral; got :$level.")
    end
end

"Convenience wrapper: `decompose(m; level=:country, …)`."
decompose_country(m::ICIOModel; perspective = :world, approach = :sink) =
    decompose(m; level = :country, perspective = perspective, approach = approach)
"Convenience wrapper: `decompose(m; level=:sector)`."
decompose_sector(m::ICIOModel) = decompose(m; level = :sector)
"Convenience wrapper: `decompose(m; level=:bilateral)`."
decompose_bilateral(m::ICIOModel) = decompose(m; level = :bilateral)
