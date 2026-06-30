# World-level perspective foreign value added (Borin-Mancini 2019 §5.1), country level.
#
# The world perspective records a unit of foreign value added only once in *total world
# exports* (rather than once in the exporting country's exports, as the exporter perspective
# does). It is available at the country level only (the icio restriction). The domestic side
# (GEXP, DC, DVA, VAX, REF, DDC) and the foreign content FC are perspective-invariant and taken
# from the source/exporter decomposition; only the FVA/FDC split changes.
#
#  * approach = :sink   — FVA recorded the *last* time the item leaves a non-origin border
#                         (Koopman-Wang-Wei 2014 logic), BM2019 eq. 54.
#  * approach = :source — FVA recorded the *first* time it is re-exported by a non-origin
#                         country (Miroudot-Ye 2018 logic), BM2019 eq. 52.

# ---------------------------------------------------------------------------
# World / sink foreign value added (BM2019 eq. 54), country level.
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

# ---------------------------------------------------------------------------
# World / source foreign value added (BM2019 eq. 52), country level.
#
# Collapsing the downstream bracket of eq. 52 with Σ_k Σ_l B_jk Y_kl = X_j and
# Σ_j Y_rj + Σ_{j≠r} A_rj X_j = (I − A_rr) X_r gives the bracket = E_sr exactly, so
#   FVA_s = Σ_{r≠s} Σ_{t≠s} V_t L_tt A_ts L_ss E_sr = Σ_{t≠s} V_t L_tt A_ts (L_ss E_s),
# i.e. foreign VA counted the first time it enters s (direct A_ts) and is re-exported by s.
# ---------------------------------------------------------------------------
function _fva_world_source(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, L, E, VLdom = m.A, m.L, m.E, m.VLdom

    fva = zeros(G)
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        etil = @view(L[srng, srng]) * @view(E[srng])   # L_ss E_s  (N-vector)
        w = @view(A[:, srng]) * etil                   # GN-vector: direct input requirement for etil
        acc = 0.0
        for a in 1:GN
            ctry(a, N) == s && continue                # keep only foreign origins t ≠ s
            acc += VLdom[a] * w[a]
        end
        fva[s] = acc
    end
    return fva
end

# World country decomposition (9 terms). Domestic side + FC from the source decomposition;
# only FVA/FDC use the chosen world approach.
function _world_country(m::ICIOModel, approach::Symbol)
    sc = _source_country(m)
    fva = approach === :sink ? _fva_world_sink(m) : _fva_world_source(m)
    fdc = sc.fc .- fva
    return (gexp = sc.gexp, dc = sc.dc, dva = sc.dva, vax = sc.vax, ref = sc.ref,
            ddc = sc.ddc, fc = sc.fc, fva = fva, fdc = fdc)
end
