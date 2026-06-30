# Sink-based country-perspective decomposition (Borin-Mancini 2019 §3.2), sector & bilateral.
#
# Same exporter (country-level) perimeter as the source decomposition, but value added is
# recorded the *last* time it leaves country s's border (double counted in prior shipments),
# rather than the first. Domestic content DC and foreign content FC are identical to the source
# breakdown at every bilateral cell (BM2019 §3.2); only the value-added vs double-counted split
# (DVA/DDC, FVA/FDC) and the VAX/REF/VAXIM components change.
#
# Key identities used (eq. 28–39):
#   * I + M_s = (I − A_ss) B_ss,  with M_s = Σ_{j≠s} A_sj B_js  and  B^{∤s}_{·s} = B_{·s}(I+M_s)^{-1}.
#   * "Ultimate" upstream output X_j^{(∤s→Y*)} = X_j − B^{∤s}_{js} E_{s*}  (eq. 32), so for a
#     demand vector d the modified inverse acts as  B^{∤s} d = B d − B_{·s}(I+M_s)^{-1}((ABd)_s − A_ss(Bd)_s).
#   * Sink bracket  Φ_sr = Y_sr + A_sr L_rr [ Σ_j Y_rj + Σ_{j≠r} A_rj X_j^{(∤s→Y*)} ];
#     DVA_sr = V_sB_ss·Φ_sr (coef VBdom), FVA_sr = Σ_{t≠s}V_tB_ts·Φ_sr (coef VBfor).
#   * REF uses the absorbed-in-s output  B^{∤s} Y_{*s}; VAX = DVA − REF.
#   * VAXIM (bilateral) uses absorbed-in-r output  Σ_{k≠s} B^{∤s}_{jk} Y_{kr}.

# Per-exporter scratch shared by both sink engines: returns closures over r giving the
# (s-sector) brackets Φ_sr (full) and Φ^ref_sr (absorbed in s).
struct _SinkExporter
    s::Int
    srng::UnitRange{Int}
    Xtil::Vector{Float64}      # X_j^{(∤s→Y*)}                 (GN)
    AXtil::Vector{Float64}     # A * Xtil                       (GN)
    Xtil_s::Vector{Float64}    # B^{∤s} Y_{*s}  (absorbed in s) (GN)
    AXtil_s::Vector{Float64}   # A * Xtil_s                     (GN)
end

function _sink_exporter(m::ICIOModel, s::Int, BFD, ABFD, Fs, Bsrng, Ass)
    A, B, X, E, FD = m.A, m.B, m.X, m.E, m.FD
    N = m.N
    srng = blockrange(s, N)
    ehat = Fs \ collect(@view E[srng])
    Xtil = X .- Bsrng * ehat
    AXtil = A * Xtil
    rhs_s = @view(ABFD[srng, s]) .- Ass * @view(BFD[srng, s])
    Xtil_s = @view(BFD[:, s]) .- Bsrng * (Fs \ collect(rhs_s))
    AXtil_s = A * Xtil_s
    return _SinkExporter(s, srng, Xtil, AXtil, Xtil_s, AXtil_s)
end

# Φ_sr (full) and Φ^ref_sr (absorbed in s), both N-vectors over s-sectors.
@inline function _sink_brackets(m::ICIOModel, e::_SinkExporter, r::Int, Yrow)
    A, L, FD, N = m.A, m.L, m.FD, m.N
    srng, rrng = e.srng, blockrange(r, N)
    Lrr = @view L[rrng, rrng]; Arr = @view A[rrng, rrng]; Asr = @view A[srng, rrng]
    Psi  = @view(Yrow[rrng]) .+ @view(e.AXtil[rrng])   .- Arr * @view(e.Xtil[rrng])
    Phi  = @view(FD[srng, r]) .+ Asr * (Lrr * Psi)
    Psir = @view(FD[rrng, e.s]) .+ @view(e.AXtil_s[rrng]) .- Arr * @view(e.Xtil_s[rrng])
    Phir = Asr * (Lrr * Psir)
    return Phi, Phir
end

# ---------------------------------------------------------------------------
# Exporter / sink, sector level (9 terms per country-sector)
# ---------------------------------------------------------------------------
function _sink_sector(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, B, E = m.A, m.B, m.E
    VBdom, VBfor = m.VBdom, m.VBfor
    BFD = _BFD(m); ABFD = A * BFD
    Yrow = vec(sum(m.FD; dims = 2))
    II_N = Matrix{Float64}(I, N, N)

    dva = zeros(GN); fva = zeros(GN); ref = zeros(GN)
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        Ass = @view A[srng, srng]; Bsrng = @view B[:, srng]
        Fs = lu((II_N .- Ass) * @view(B[srng, srng]))
        e = _sink_exporter(m, s, BFD, ABFD, Fs, Bsrng, Ass)
        for r in 1:G
            r == s && continue
            Phi, Phir = _sink_brackets(m, e, r, Yrow)
            for (k, n) in enumerate(srng)
                dva[n] += VBdom[n] * Phi[k]
                fva[n] += VBfor[n] * Phi[k]
                ref[n] += VBdom[n] * Phir[k]
            end
        end
    end
    gexp = copy(E); dc = VBdom .* E; fc = VBfor .* E
    ddc = dc .- dva; fdc = fc .- fva; vax = dva .- ref
    return (gexp = gexp, dc = dc, dva = dva, vax = vax, ref = ref,
            ddc = ddc, fc = fc, fva = fva, fdc = fdc)
end

# ---------------------------------------------------------------------------
# Exporter / sink, bilateral-sector level (10 terms per s,n,r with r≠s; adds VAXIM)
# ---------------------------------------------------------------------------
function _sink_bilateral(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, B, L, ESR, FD = m.A, m.B, m.L, m.ESR, m.FD
    VBdom, VBfor = m.VBdom, m.VBfor
    BFD = _BFD(m); ABFD = A * BFD
    Yrow = vec(sum(FD; dims = 2))
    II_N = Matrix{Float64}(I, N, N)

    nrow = G * (G - 1) * N
    exp_g = Vector{Int}(undef, nrow); exp_n = Vector{Int}(undef, nrow); imp_r = Vector{Int}(undef, nrow)
    gexp = Vector{Float64}(undef, nrow); dc = similar(gexp); dva = similar(gexp)
    vax = similar(gexp); vaxim = similar(gexp); ref = similar(gexp); ddc = similar(gexp)
    fc = similar(gexp); fva = similar(gexp); fdc = similar(gexp)

    row = 0
    @inbounds for s in 1:G
        srng = blockrange(s, N)
        Ass = @view A[srng, srng]; Bsrng = @view B[:, srng]
        Fs = lu((II_N .- Ass) * @view(B[srng, srng]))
        e = _sink_exporter(m, s, BFD, ABFD, Fs, Bsrng, Ass)
        for r in 1:G
            r == s && continue
            rrng = blockrange(r, N)
            Lrr = @view L[rrng, rrng]; Arr = @view A[rrng, rrng]; Asr = @view A[srng, rrng]
            Phi, Phir = _sink_brackets(m, e, r, Yrow)
            # VAXIM: absorbed-in-r output of each origin, Σ_{k≠s} B^{∤s}_{jk} Y_{kr}
            xr = @view(BFD[:, r]) .- Bsrng *
                 (Fs \ collect(@view(ABFD[srng, r]) .- Ass * @view(BFD[srng, r]) .+ @view(FD[srng, r])))
            arx = @view(A[rrng, :]) * xr .- Arr * @view(xr[rrng])      # Σ_{j≠r} A_rj xr_j (N, r-sectors)
            Psivx = @view(FD[rrng, r]) .+ arx
            Phivx = @view(FD[srng, r]) .+ Asr * (Lrr * Psivx)
            for (k, n) in enumerate(srng)
                row += 1
                e_sr = ESR[n, r]
                DCv = VBdom[n] * e_sr; FCv = VBfor[n] * e_sr
                DVAv = VBdom[n] * Phi[k]; FVAv = VBfor[n] * Phi[k]
                REFv = VBdom[n] * Phir[k]
                exp_g[row] = s; exp_n[row] = k; imp_r[row] = r
                gexp[row] = e_sr; dc[row] = DCv; fc[row] = FCv
                dva[row] = DVAv; ddc[row] = DCv - DVAv
                fva[row] = FVAv; fdc[row] = FCv - FVAv
                ref[row] = REFv; vax[row] = DVAv - REFv
                vaxim[row] = VBdom[n] * Phivx[k]
            end
        end
    end
    return (exp_g = exp_g, exp_n = exp_n, imp_r = imp_r,
            terms = (gexp = gexp, dc = dc, dva = dva, vax = vax, vaxim = vaxim, ref = ref,
                     ddc = ddc, fc = fc, fva = fva, fdc = fdc))
end
