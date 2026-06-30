# Importer-perspective decomposition of gross imports (Borin-Mancini 2019 §4.3, eq. 51).
#
# The perimeter is the importing country r's border as a whole: gross imports of r are split
# into value added (VA, GDP of some origin j embodied in r's imports, counted once at r's
# border) and double counting (DC). There is no domestic/foreign distinction (the notion is
# relative to the exporter), but VA and DC can be attributed to a country of origin j.
#
#   u_N E_{*r} = Σ_j V_j Σ_{s≠r} B̃^r_{js} E_{sr}        (VA)
#              + Σ_j V_j Σ_{t≠r} B̃^r_{jt} Σ_{s≠r} A_tr B_rs E_{sr}   (DC)
#
# where B̃^r = (I − Ã^r)^{-1} and Ã^r zeroes the off-diagonal blocks A_tr (t≠r) of the r-th
# block-column. This is a rank-N (column-block r) change of (I−A), so B̃^r follows from a
# Woodbury update of the cached B:
#   I_N + (B M_col)_{rr} = B_rr (I − A_rr) ≡ K_r,   B M_col = B_{·r}(I−A_rr) − I_{·r},
#   B̃^r d = B d − (B M_col) K_r^{-1} (B d)_r .
#
# Let d^r = ESR[:, r] (imports of r by origin country-sector). Then VA = V'·(B̃^r d^r); the
# inner re-import demand is h^r = (A_{·r} (B d^r)_r) with the r-block zeroed, and DC = V'·(B̃^r h^r).

# Apply B̃^r to a demand vector via the Woodbury update. `Bd` = B*d (GN), `BMcol` = B M_col (GN×N),
# `Fr` = lu(K_r), `rrng` = r-block.
@inline _btilde_apply(Bd, BMcol, Fr, rrng) = Bd .- BMcol * (Fr \ collect(@view Bd[rrng]))

# Core: returns per-importer, per-origin VA and DC (GN×G arrays indexed [origin-sector, importer]).
function _imports_core(m::ICIOModel)
    G, N, GN = m.G, m.N, m.GN
    A, B, Vc, ESR = m.A, m.B, m.V, m.ESR
    BESR = B * ESR                       # GN×G, BESR[:,r] = B d^r
    II_N = Matrix{Float64}(I, N, N)

    va = zeros(GN, G)   # va[a, r] = Vc[a]·(B̃^r d^r)[a]  (origin-sector a, importer r)
    dc = zeros(GN, G)
    @inbounds for r in 1:G
        rrng = blockrange(r, N)
        Irr = II_N .- @view A[rrng, rrng]
        Fr = lu(@view(B[rrng, rrng]) * Irr)              # K_r = B_rr (I − A_rr)
        BMcol = @view(B[:, rrng]) * Irr                  # B_{·r}(I − A_rr)
        @views BMcol[rrng, :] .-= II_N                   # − I_{·r}
        Bdr = collect(@view BESR[:, r])
        btil_dr = _btilde_apply(Bdr, BMcol, Fr, rrng)    # B̃^r d^r
        gr = @view Bdr[rrng]                             # g^r = (B d^r)_r
        hr = @view(A[:, rrng]) * gr                      # A_{·r} g^r
        @views hr[rrng] .= 0.0                           # keep t ≠ r
        Bhr = B * hr
        btil_hr = _btilde_apply(Bhr, BMcol, Fr, rrng)    # B̃^r h^r
        for a in 1:GN
            va[a, r] = Vc[a] * btil_dr[a]
            dc[a, r] = Vc[a] * btil_hr[a]
        end
    end
    return va, dc
end

# Country level: one row per importer r (gimp, va, dc).
function _imports_country(m::ICIOModel)
    G, N = m.G, m.N
    va, dc = _imports_core(m)
    VA = vec(sum(va; dims = 1)); DC = vec(sum(dc; dims = 1))
    return (gimp = VA .+ DC, va = VA, dc = DC)
end

# Bilateral level: VA-origin breakdown of r's imports — one row per (importer r, origin j),
# all j (origin j = r is the re-imported domestic value added). va/dc sum over j to total
# gross imports of r, so this is a complete partition of imports by country of origin.
function _imports_bilateral(m::ICIOModel)
    G, N = m.G, m.N
    va, dc = _imports_core(m)
    nrow = G * G
    imp_r = Vector{Int}(undef, nrow); ori_j = Vector{Int}(undef, nrow)
    VA = Vector{Float64}(undef, nrow); DC = similar(VA)
    row = 0
    @inbounds for r in 1:G
        for j in 1:G
            row += 1
            jrng = blockrange(j, N)
            imp_r[row] = r; ori_j[row] = j
            v = 0.0; d = 0.0
            for a in jrng
                v += va[a, r]; d += dc[a, r]
            end
            VA[row] = v; DC[row] = d
        end
    end
    return (imp_r = imp_r, ori_j = ori_j, terms = (va = VA, dc = DC))
end
