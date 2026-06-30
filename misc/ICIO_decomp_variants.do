*** Stata `icio` reference outputs for the extended GlobalValueChains.jl decomposition variants ***
* Companion to misc/ICIO_decomp_variants.jl. Produces *_STATA.csv reference files for the new
* perspectives/approaches so the Julia harness can diff against them (expect ~1e-7 agreement).
*
* Variants produced (one year, 2015):
*   country   world / source              -> EM_GVC_KWW_WS_BM19_STATA.csv     (9 terms)
*   sector    exporter / sink             -> EM_GVC_SEC_SINK_BM19_STATA.csv   (9 terms)
*   bilateral exporter / sink (sample)    -> EM_GVC_BIL_SINK_SAMPLE_STATA.csv (10 terms, incl VAXIM)
*   imports   importer perspective        -> EM_GVC_IMP_BM19_STATA.csv        (gimp va dc)
*
* Run from the EMERGING project root (parent of ICIO_CSV/):
*   do path/to/GlobalValueChains.jl/misc/ICIO_decomp_variants.do
* Then re-run the Julia harness to diff:
*   julia --project=/path/to/GlobalValueChains.jl misc/ICIO_decomp_variants.jl

clear
global data EM
cap global datadir "/Users/sebastiankrantz/Documents/Data/EMERGING/ICIO_CSV/EMERGING_Broad_Sectors"
if "$datadir" == "" global datadir "ICIO_CSV/EMERGING_Broad_Sectors"
global csv_path "$datadir"
global year 2015

* sector count
import delimited "${csv_path}/${data}_countrylist.csv", delimiters(",") varnames(nonames)
global nctry = _N
levelsof v1, local(allctry) clean
global allctry "`allctry'"      // promote to global so it survives the `clear`s below
import delimited "${csv_path}/${data}_${year}.csv", delimiters(",") varnames(nonames) clear
global nsec = _N / $nctry
clear

* bilateral sample (keep in sync with EXPORTERS/IMPORTERS in the .jl)
global exporters CHN DEU USA ZAF IND
global importers USA CHN DEU NGA

icio_clean
icio_load, iciot(user, userp($csv_path) tablen(${data}_${year}.csv) countrylist(${data}_countrylist.csv))
qui icio, info

* ---------------------------------------------------------------------------
* 1) Country-level, world / source  (9 terms: GEXP DC DVA VAX REF DDC FC FVA FDC)
* ---------------------------------------------------------------------------
clear
set obs $nctry
gen country = ""
foreach v in gexp dc dva vax ref ddc fc fva fdc {
    gen `v' = .
}
local i = 0
foreach ce of global allctry {
    local i = `i' + 1
    qui icio, exporter(`ce') perspective(world) approach(source) output(detailed)
    mat res = r(detailed)
    qui replace country = "`ce'" in `i'
    qui replace gexp = res[1,1] in `i'
    qui replace dc   = res[2,1] in `i'
    qui replace dva  = res[3,1] in `i'
    qui replace vax  = res[4,1] in `i'
    qui replace ref  = res[5,1] in `i'
    qui replace ddc  = res[6,1] in `i'
    qui replace fc   = res[7,1] in `i'
    qui replace fva  = res[8,1] in `i'
    qui replace fdc  = res[9,1] in `i'
}
outsheet using "${csv_path}/${data}_GVC_KWW_WS_BM19_STATA.csv", comma replace
di "Saved EM_GVC_KWW_WS_BM19_STATA.csv"

* ---------------------------------------------------------------------------
* 2) Sector-level, exporter / sink  (9 terms; one block of $nsec columns per exporter)
* ---------------------------------------------------------------------------
clear
set obs `=$nctry * $nsec'
gen from_region = ""
gen from_sector = .
foreach v in gexp dc dva vax ref ddc fc fva fdc {
    gen `v' = .
}
local blk = 0
foreach ce of global allctry {
    qui icio, exporter(`ce', all) perspective(exporter) approach(sink) output(detailed)
    mat res = r(detailed)                       // 9 rows x nsec cols
    forvalues sct = 1/$nsec {
        local row = `blk' * $nsec + `sct'
        qui replace from_region = "`ce'" in `row'
        qui replace from_sector = `sct' in `row'
        qui replace gexp = res[1,`sct'] in `row'
        qui replace dc   = res[2,`sct'] in `row'
        qui replace dva  = res[3,`sct'] in `row'
        qui replace vax  = res[4,`sct'] in `row'
        qui replace ref  = res[5,`sct'] in `row'
        qui replace ddc  = res[6,`sct'] in `row'
        qui replace fc   = res[7,`sct'] in `row'
        qui replace fva  = res[8,`sct'] in `row'
        qui replace fdc  = res[9,`sct'] in `row'
    }
    local blk = `blk' + 1
}
outsheet using "${csv_path}/${data}_GVC_SEC_SINK_BM19_STATA.csv", comma replace
di "Saved EM_GVC_SEC_SINK_BM19_STATA.csv"

* ---------------------------------------------------------------------------
* 3) Bilateral sector-level, exporter / sink (10 terms incl VAXIM), sample pairs
*    row order: GEXP DC DVA VAX VAXIM REF DDC FC FVA FDC
* ---------------------------------------------------------------------------
clear
local npairs = 0
foreach ci in $importers {
    foreach ce in $exporters {
        if "`ce'" != "`ci'" local npairs = `npairs' + 1
    }
}
set obs `=`npairs' * $nsec'
gen from_region = ""
gen from_sector = .
gen to_region = ""
foreach v in gexp dc dva vax vaxim ref ddc fc fva fdc {
    gen `v' = .
}
local iter = 0
foreach ci in $importers {
    foreach ce in $exporters {
        if "`ce'" == "`ci'" continue
        qui icio, exporter(`ce', all) importer(`ci') perspective(exporter) approach(sink) output(detailed)
        mat res = r(detailed)
        local iter = `iter' + 1
        forvalues sct = 1/$nsec {
            local row = (`iter'-1)*$nsec + `sct'
            qui replace from_region = "`ce'" in `row'
            qui replace from_sector = `sct' in `row'
            qui replace to_region = "`ci'" in `row'
            qui replace gexp  = res[1,`sct'] in `row'
            qui replace dc    = res[2,`sct'] in `row'
            qui replace dva   = res[3,`sct'] in `row'
            qui replace vax   = res[4,`sct'] in `row'
            qui replace vaxim = res[5,`sct'] in `row'
            qui replace ref   = res[6,`sct'] in `row'
            qui replace ddc   = res[7,`sct'] in `row'
            qui replace fc    = res[8,`sct'] in `row'
            qui replace fva   = res[9,`sct'] in `row'
            qui replace fdc   = res[10,`sct'] in `row'
        }
    }
}
outsheet using "${csv_path}/${data}_GVC_BIL_SINK_SAMPLE_STATA.csv", comma replace
di "Saved EM_GVC_BIL_SINK_SAMPLE_STATA.csv"

* ---------------------------------------------------------------------------
* 4) Imports, importer perspective (gross imports split into VA + DC), by country
*    NOTE: icio returns gross imports and total VA; DC = gimp - va.
* ---------------------------------------------------------------------------
clear
set obs $nctry
gen importer = ""
gen gimp = .
gen va = .
gen dc = .
local i = 0
foreach ci of global allctry {
    local i = `i' + 1
    qui icio, importer(`ci') perspective(importer) output(gtrade)
    mat g = r(gtrade)
    qui icio, importer(`ci') perspective(importer) output(va)
    mat v = r(va)
    qui replace importer = "`ci'" in `i'
    qui replace gimp = g[1,1] in `i'
    qui replace va   = v[1,1] in `i'
    qui replace dc   = g[1,1] - v[1,1] in `i'
}
outsheet using "${csv_path}/${data}_GVC_IMP_BM19_STATA.csv", comma replace
di "Saved EM_GVC_IMP_BM19_STATA.csv"

* ---------------------------------------------------------------------------
* 5) Sector-level, self (sectexp) perimeter  (9 terms: GEXP DC DVA VAX REF DDC FC FVA FDC)
* ---------------------------------------------------------------------------
clear
set obs `=$nctry * $nsec'
gen from_region = ""
gen from_sector = .
foreach v in gexp dc dva vax ref ddc fc fva fdc {
    gen `v' = .
}
local blk = 0
foreach ce of global allctry {
    qui icio, exporter(`ce', all) perspective(sectexp) output(detailed)
    mat res = r(detailed)                       // 9 rows x nsec cols
    forvalues sct = 1/$nsec {
        local row = `blk' * $nsec + `sct'
        qui replace from_region = "`ce'" in `row'
        qui replace from_sector = `sct' in `row'
        qui replace gexp = res[1,`sct'] in `row'
        qui replace dc   = res[2,`sct'] in `row'
        qui replace dva  = res[3,`sct'] in `row'
        qui replace vax  = res[4,`sct'] in `row'
        qui replace ref  = res[5,`sct'] in `row'
        qui replace ddc  = res[6,`sct'] in `row'
        qui replace fc   = res[7,`sct'] in `row'
        qui replace fva  = res[8,`sct'] in `row'
        qui replace fdc  = res[9,`sct'] in `row'
    }
    local blk = `blk' + 1
}
outsheet using "${csv_path}/${data}_GVC_SEC_SELF_BM19_STATA.csv", comma replace
di "Saved EM_GVC_SEC_SELF_BM19_STATA.csv"

* ---------------------------------------------------------------------------
* 6) Bilateral sector-level, self (sectbil) perimeter (9 terms), sample pairs
*    row order: GEXP DC DVA VAX REF DDC FC FVA FDC
* ---------------------------------------------------------------------------
clear
local npairs = 0
foreach ci in $importers {
    foreach ce in $exporters {
        if "`ce'" != "`ci'" local npairs = `npairs' + 1
    }
}
set obs `=`npairs' * $nsec'
gen from_region = ""
gen from_sector = .
gen to_region = ""
foreach v in gexp dc dva vax ref ddc fc fva fdc {
    gen `v' = .
}
local iter = 0
foreach ci in $importers {
    foreach ce in $exporters {
        if "`ce'" == "`ci'" continue
        qui icio, exporter(`ce', all) importer(`ci') perspective(sectbil) output(detailed)
        mat res = r(detailed)
        local iter = `iter' + 1
        forvalues sct = 1/$nsec {
            local row = (`iter'-1)*$nsec + `sct'
            qui replace from_region = "`ce'" in `row'
            qui replace from_sector = `sct' in `row'
            qui replace to_region = "`ci'" in `row'
            qui replace gexp = res[1,`sct'] in `row'
            qui replace dc   = res[2,`sct'] in `row'
            qui replace dva  = res[3,`sct'] in `row'
            qui replace vax  = res[4,`sct'] in `row'
            qui replace ref  = res[5,`sct'] in `row'
            qui replace ddc  = res[6,`sct'] in `row'
            qui replace fc   = res[7,`sct'] in `row'
            qui replace fva  = res[8,`sct'] in `row'
            qui replace fdc  = res[9,`sct'] in `row'
        }
    }
}
outsheet using "${csv_path}/${data}_GVC_BIL_SELF_SAMPLE_STATA.csv", comma replace
di "Saved EM_GVC_BIL_SELF_SAMPLE_STATA.csv"

di "Done. Re-run misc/ICIO_decomp_variants.jl to diff."
