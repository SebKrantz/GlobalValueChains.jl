# Self-perimeter decompositions (Borin-Mancini 2019 §4): the perimeter for "double counting"
# is the export flow itself, not the exporting country as a whole. This yields the broader
# domestic value-added notion of Johnson (2018) / Los et al. (2016): DVA★ ⊇ DVAsource, DVAsink.
#
#  * level = :sector    (sectexp, §4.3) — perimeter = a sector's total exports; zero a_{sj,n} ∀ j≠s.
#  * level = :bilateral (sectbil, eq. 47–49) — perimeter = one sectoral-bilateral flow; zero a_{sr,n·}.
#
# Both modifications change a single row (s,n) of A, so the modified Leontief inverse follows
# from a rank-1 Woodbury update of the cached B:  B^{∗}_{·,(s,n)} = B_{·,(s,n)} / (1 + α),
# α = P[n,n], with the "re-import" coefficient P = A_sr B_rs (sectbil) or P = M_s = Σ_{j≠s} A_sj B_js
# (sectexp). The self-perimeter value-added terms (eq. 49 terms 1 and 3) are then
#   DVA★ = V_s B^{∗}_ss E = VBdom_{(s,n)}/(1+α) · e,   FVA★ = Σ_{t≠s}V_t B^{∗}_ts E = VBfor_{(s,n)}/(1+α) · e.
# Domestic and foreign *content* are perimeter-invariant (BM2019 §3.2), so DC = VBdom·e,
# FC = VBfor·e, and the double-counted terms close the decomposition: DDC = DC − DVA★,
# FDC = FC − FVA★ (guaranteeing GEXP = DC + FC = DVA★+DDC+FVA★+FDC). VAX/REF (the abroad/home
# split of DVA★) are perimeter-specific and only provided under perspective = :exporter.

# (DVA★, FVA★) for one export cell. `gi` = global index of (s,n); `den = 1 + P[n,n]`.
@inline _self_dva_fva(VBdom, VBfor, gi::Int, den::Float64, e::Float64) =
    (VBdom[gi] / den * e, VBfor[gi] / den * e)

# ---------------------------------------------------------------------------
# Sectoral-exporter (sectexp) self perimeter — one row per country-sector (7 terms)
# ---------------------------------------------------------------------------
function _self_sector(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, B, E = m.A, m.B, m.E
    VBdom, VBfor = m.VBdom, m.VBfor
    II_N = Matrix{Float64}(I, N, N)

    dva = zeros(GN); fva = zeros(GN)
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        Bss = @view B[srng, srng]
        Ms = (Bss .- II_N) .- @view(A[srng, srng]) * Bss        # M_s = Σ_{j≠s} A_sj B_js
        for n in 1:N
            gi = (s - 1) * N + n
            dva[gi], fva[gi] = _self_dva_fva(VBdom, VBfor, gi, 1.0 + Ms[n, n], E[gi])
        end
    end
    gexp = copy(E); dc = VBdom .* E; fc = VBfor .* E
    ddc = dc .- dva; fdc = fc .- fva
    return (gexp = gexp, dc = dc, dva = dva, ddc = ddc, fc = fc, fva = fva, fdc = fdc)
end

# ---------------------------------------------------------------------------
# Sectoral-bilateral (sectbil) self perimeter — one row per (s,n,r), r≠s (7 terms)
# ---------------------------------------------------------------------------
function _self_bilateral(m::ICIOModel)
    G, N = m.G, m.N
    A, B, ESR = m.A, m.B, m.ESR
    VBdom, VBfor = m.VBdom, m.VBfor

    nrow = G * (G - 1) * N
    exp_g = Vector{Int}(undef, nrow); exp_n = Vector{Int}(undef, nrow); imp_r = Vector{Int}(undef, nrow)
    gexp = Vector{Float64}(undef, nrow); dc = similar(gexp); dva = similar(gexp); ddc = similar(gexp)
    fc = similar(gexp); fva = similar(gexp); fdc = similar(gexp)

    row = 0
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        for r in 1:G
            r == s && continue
            rrng = blockrange(r, N)
            Psr = @view(A[srng, rrng]) * @view(B[rrng, srng])     # A_sr B_rs  (N×N)
            for n in 1:N
                row += 1
                gi = (s - 1) * N + n
                e = ESR[gi, r]
                d, f = _self_dva_fva(VBdom, VBfor, gi, 1.0 + Psr[n, n], e)
                DCv = VBdom[gi] * e; FCv = VBfor[gi] * e
                exp_g[row] = s; exp_n[row] = n; imp_r[row] = r
                gexp[row] = e; dc[row] = DCv; fc[row] = FCv
                dva[row] = d; ddc[row] = DCv - d; fva[row] = f; fdc[row] = FCv - f
            end
        end
    end
    return (exp_g = exp_g, exp_n = exp_n, imp_r = imp_r,
            terms = (gexp = gexp, dc = dc, dva = dva, ddc = ddc, fc = fc, fva = fva, fdc = fdc))
end
