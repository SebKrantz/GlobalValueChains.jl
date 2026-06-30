*** Bilateral Sector-Level Decomposition — VALIDATION SAMPLE ***
* Mirrors the third block of ICIO_decomp.do, but restricted to a handful of
* exporter/importer pairs (and one year) so it runs in seconds. Produces a
* reference file to compare against GlobalValueChains.jl's `decompose(m; level = :bilateral)`.
*
* Data layout (same as misc/ICIO_decomp_bil_sample.jl):
*   cd to the EMERGING project root (parent of ICIO_CSV/), or set global datadir below.
*   Julia equivalent:  export ICIO_DATA=/path/to/EMERGING_Broad_Sectors
*
* From the EMERGING project root:  do path/to/GlobalValueChains.jl/misc/ICIO_decomp_bil_sample.do
* Then compare:  julia --project=/path/to/GlobalValueChains.jl misc/ICIO_decomp_bil_sample.jl

clear
global data EM
cap global datadir ""
if "$datadir" == "" global datadir "ICIO_CSV/EMERGING_Broad_Sectors"
global csv_path "$datadir"

* Country count and sector count (as in ICIO_decomp.do SETUP)
import delimited "${csv_path}/${data}_countrylist.csv", delimiters(",") varnames(nonames)
global nctry = _N
import delimited "${csv_path}/${data}_2015.csv", delimiters(",") varnames(nonames) clear
global nsec = _N / $nctry
clear

di "Sectors: ${nsec}"

* ---- Validation sample: a few diverse exporters and importers, one year ----
global years 2015   // keep in sync with SAMPLE_YEAR in misc/ICIO_decomp_bil_sample.jl
global exporters CHN DEU USA ZAF ABW IND
global importers USA CHN DEU NGA

* Count valid (exporter != importer) pairs to size the dataset
local npairs = 0
foreach ci in $importers {
	foreach ce in $exporters {
		if "`ce'" != "`ci'" local npairs = `npairs' + 1
	}
}
local nobs = `npairs' * $nsec
set obs `nobs'

* ID + term variables (same 13 terms / row order as block 3 of ICIO_decomp.do)
gen year = .
gen from_region = ""
gen from_sector = .
gen to_region = ""
foreach var in gexp dc dva vax davax ref ddc fc fva fdc gvc gvcb gvcf {
	gen `var' = .
}

local iter = 0
foreach y in $years {

	di "Year: `y'"
	icio_clean
	icio_load, iciot(user, userp($csv_path) tablen(${data}_`y'.csv) countrylist(${data}_countrylist.csv))
	qui icio, info

	foreach ci in $importers {
		foreach ce in $exporters {
			if "`ce'" == "`ci'" continue
			qui icio, exporter(`ce', all) importer(`ci') perspective(exporter) approach(source) output(detailed)
			mat res = r(detailed)
			local iter = `iter' + 1
			local start = (`iter'-1)*$nsec + 1
			local end = `start' + $nsec - 1
			local i = 0
			qui forvalues j = `start'/`end' {
				local i = `i' + 1
				replace year = `y' in `j'
				replace from_region = "`ce'" in `j'
				replace from_sector = `i' in `j'
				replace to_region = "`ci'" in `j'
				replace gexp  = res[1,  `i'] in `j'
				replace dc    = res[2,  `i'] in `j'
				replace dva   = res[3,  `i'] in `j'
				replace vax   = res[4,  `i'] in `j'
				replace davax = res[5,  `i'] in `j'
				replace ref   = res[6,  `i'] in `j'
				replace ddc   = res[7,  `i'] in `j'
				replace fc    = res[8,  `i'] in `j'
				replace fva   = res[9,  `i'] in `j'
				replace fdc   = res[10, `i'] in `j'
				replace gvc   = res[11, `i'] in `j'
				replace gvcb  = res[12, `i'] in `j'
				replace gvcf  = res[13, `i'] in `j'
			}
		}
	}
}

outsheet using "${csv_path}/${data}_GVC_BIL_SEC_SAMPLE.csv", comma replace
di "Saved ${csv_path}/${data}_GVC_BIL_SEC_SAMPLE.csv"
