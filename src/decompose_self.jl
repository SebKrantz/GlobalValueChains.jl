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
# (sectexp). The self value-added terms (eq. 49 terms 1 and 3) are then DVA★ = VBdom_{(s,n)}/(1+α)·e
# and FVA★ = VBfor_{(s,n)}/(1+α)·e. Domestic and foreign *content* are perimeter-invariant
# (BM2019 §3.2), so DC = VBdom·e, FC = VBfor·e, and DDC = DC−DVA★, FDC = FC−FVA★ close the
# decomposition (GEXP = DC + FC).
#
# VAX/REF split DVA★ by where the export is ultimately absorbed (abroad vs. back home in s). This
# is a *downstream* property of the flow — every unit of value added in the export travels the
# same forward path regardless of the upstream perimeter — so the reflection share is
# perimeter-invariant: VAX★ = (DVA★/e)·VAXE, where VAXE is the export's downstream-absorbed-abroad
# demand (the same expansion used by the exporter/source decomposition, eq. 19/20). Hence
# VAX★/DVA★ = VAXsource/DVAsource and REF★ = DVA★ − VAX★.

# ---------------------------------------------------------------------------
# Sectoral-exporter (sectexp) self perimeter — one row per country-sector (9 terms)
# ---------------------------------------------------------------------------
function _self_sector(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, B, X, E, FD = m.A, m.B, m.X, m.E, m.FD
    VBdom, VBfor = m.VBdom, m.VBfor
    BFD = _BFD(m)
    II_N = Matrix{Float64}(I, N, N)

    VAXEsum = zeros(GN)           # Σ_{r≠s} (Y_sr + A_sr (X_r − BFD[r,s]))  (export absorbed abroad)
    dva = zeros(GN); vax = zeros(GN); fva = zeros(GN)
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        for r in 1:G
            r == s && continue
            rrng = blockrange(r, N)
            vaxe = @view(A[srng, rrng]) * (@view(X[rrng]) .- @view(BFD[rrng, s]))   # N-vector
            for (k, gi) in enumerate(srng)
                VAXEsum[gi] += FD[gi, r] + vaxe[k]
            end
        end
        Bss = @view B[srng, srng]
        Ms = (Bss .- II_N) .- @view(A[srng, srng]) * Bss        # M_s = Σ_{j≠s} A_sj B_js
        for n in 1:N
            gi = (s - 1) * N + n
            cstar = VBdom[gi] / (1.0 + Ms[n, n])                # self DVA coefficient
            dva[gi] = cstar * E[gi]
            vax[gi] = cstar * VAXEsum[gi]
            fva[gi] = VBfor[gi] / (1.0 + Ms[n, n]) * E[gi]
        end
    end
    gexp = copy(E); dc = VBdom .* E; fc = VBfor .* E
    ref = dva .- vax; ddc = dc .- dva; fdc = fc .- fva
    return (gexp = gexp, dc = dc, dva = dva, vax = vax, ref = ref,
            ddc = ddc, fc = fc, fva = fva, fdc = fdc)
end

# ---------------------------------------------------------------------------
# Sectoral-bilateral (sectbil) self perimeter — one row per (s,n,r), r≠s (9 terms)
# ---------------------------------------------------------------------------
function _self_bilateral(m::ICIOModel)
    G, N = m.G, m.N
    A, B, X, ESR, FD = m.A, m.B, m.X, m.ESR, m.FD
    VBdom, VBfor = m.VBdom, m.VBfor
    BFD = _BFD(m)

    nrow = G * (G - 1) * N
    exp_g = Vector{Int}(undef, nrow); exp_n = Vector{Int}(undef, nrow); imp_r = Vector{Int}(undef, nrow)
    gexp = Vector{Float64}(undef, nrow); dc = similar(gexp); dva = similar(gexp)
    vax = similar(gexp); ref = similar(gexp); ddc = similar(gexp)
    fc = similar(gexp); fva = similar(gexp); fdc = similar(gexp)

    row = 0
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        for r in 1:G
            r == s && continue
            rrng = blockrange(r, N)
            Psr = @view(A[srng, rrng]) * @view(B[rrng, srng])                       # A_sr B_rs  (N×N)
            vaxe = @view(A[srng, rrng]) * (@view(X[rrng]) .- @view(BFD[rrng, s]))   # N-vector
            for n in 1:N
                row += 1
                gi = (s - 1) * N + n
                e = ESR[gi, r]
                cstar = VBdom[gi] / (1.0 + Psr[n, n])
                DVAv = cstar * e
                VAXv = cstar * (FD[gi, r] + vaxe[n])           # downstream-abroad share of e
                DCv = VBdom[gi] * e; FCv = VBfor[gi] * e
                FVAv = VBfor[gi] / (1.0 + Psr[n, n]) * e
                exp_g[row] = s; exp_n[row] = n; imp_r[row] = r
                gexp[row] = e; dc[row] = DCv; fc[row] = FCv
                dva[row] = DVAv; vax[row] = VAXv; ref[row] = DVAv - VAXv
                ddc[row] = DCv - DVAv; fva[row] = FVAv; fdc[row] = FCv - FVAv
            end
        end
    end
    return (exp_g = exp_g, exp_n = exp_n, imp_r = imp_r,
            terms = (gexp = gexp, dc = dc, dva = dva, vax = vax, ref = ref,
                     ddc = ddc, fc = fc, fva = fva, fdc = fdc))
end
