*** SETUP: Define Global Variables ***
*
* Data layout (same as misc/ICIO_decomp.jl):
*   cd to the EMERGING project root (parent of ICIO_CSV/), or set global datadir below.
*   Julia equivalent:  export ICIO_DATA=/path/to/EMERGING_Broad_Sectors
*
* From the EMERGING project root:  do path/to/ICIO.jl/misc/ICIO_decomp.do
* Or set global datadir before running (see Julia's ICIO_DATA).

clear
global data EM
cap global datadir ""
if "$datadir" == "" global datadir "ICIO_CSV/EMERGING_Broad_Sectors"
global csv_path "$datadir"
import delimited "${csv_path}/${data}_countrylist.csv", delimiters(",") varnames(nonames)
levelsof v1, clean
global countries `r(levels)'
global nctry = _N
import delimited "${csv_path}/${data}_2015.csv", delimiters(",")  varnames(nonames) clear
global nsec = _N / $nctry
clear
global years 2015 2018 2021 2023   // keep in sync with YEARS in misc/ICIO_decomp.jl
global nyears 4

di "Dataset: ${data}, Countries: ${nctry}, Sectors: ${nsec}, Years: ${nyears}"


*** KWW Decomposition (Corrected) ***

clear 
local nobs = $nyears * $nctry
set obs `nobs'

* ID Variables
gen year = .
gen country = ""

* Generate Additional Names
foreach var in gexp dc dva vax ref ddc fc fva fdc {
	gen `var' = .
}

foreach y in $years {
	
	di "Year: `y'"
	
	* Load ICIO Table
	icio_clean
	icio_load, iciot(user, userp($csv_path) tablen(${data}_`y'.csv) countrylist(${data}_countrylist.csv))  
	icio, info
	
	qui foreach c in $countries {
		icio, exporter(`c') perspective(world) approach(sink)
		// matlist r(detailed)
		mat res = r(detailed)
		local j = `j'+1
		replace year = `y' in `j'
		replace country = "`c'" in `j'
		replace gexp = res[1, 1] in `j'
		replace dc = res[2, 1] in `j'
		replace dva = res[3, 1] in `j'
		replace vax = res[4, 1] in `j'
		replace ref = res[5, 1] in `j'
		replace ddc = res[6, 1] in `j'
		replace fc = res[7, 1] in `j'
		replace fva = res[8, 1] in `j'
		replace fdc = res[9, 1] in `j'
	}
}

outsheet using "${csv_path}/${data}_GVC_KWW_BM19.csv", comma replace




*** Exporter-Sector Level Decomposition ***

clear 
local nobs = $nyears * $nctry * $nsec
set obs `nobs'

* ID Variables
gen year = .
gen from_region = ""
gen from_sector = .

* Generate Additional Names
foreach var in gexp dc dva vax davax ref ddc fc fva fdc gvc gvcb gvcf {
	gen `var' = .
}
	
foreach y in $years {
	
	di "Year: `y'"
	
	* Load ICIO Table
	icio_clean
	icio_load, iciot(user, userp($csv_path) tablen(${data}_`y'.csv) countrylist(${data}_countrylist.csv))  
	icio, info

	* Detailed Decomposition
	foreach ce in $countries {
		qui icio, exporter(`ce', all) perspective(exporter) approach(source) output(detailed)
		// matlist r(detailed)
		mat res = r(detailed)
		local iter = `iter'+1
		local start = (`iter'-1)*$nsec + 1
		local end = `start' + $nsec - 1
		* di "EXP: `ce', IMP: `ci'-`si', Rows: `start'-`end'"
		local i = 0
		qui forvalues j = `start'/`end' {
			local i = `i'+1
			replace year = `y' in `j'
			replace from_region = "`ce'" in `j'
			replace from_sector = `i' in `j'
			replace gexp = res[1, `i'] in `j'
			replace dc = res[2, `i'] in `j'
			replace dva = res[3, `i'] in `j'
			replace vax = res[4, `i'] in `j'
			replace davax = res[5, `i'] in `j'
			replace ref = res[6, `i'] in `j'
			replace ddc = res[7, `i'] in `j'
			replace fc = res[8, `i'] in `j'
			replace fva = res[9, `i'] in `j'
			replace fdc = res[10, `i'] in `j'
			replace gvc = res[11, `i'] in `j'
			replace gvcb = res[12, `i'] in `j'
			replace gvcf = res[13, `i'] in `j'
		}			
	}
}


outsheet using "${csv_path}/${data}_GVC_SEC_BM19.csv", comma replace



*** Detailed Bilateral-Sector Decomposition ***
/* Note: Cannot have importing country sector because final demand is not available at the sector level 
   Also Note: Computations can take extremely long
*/

clear 
local nobs = $nyears * $nctry * $nsec * ($nctry-1)
set obs `nobs'

* ID Variables
gen year = .
gen from_region = ""
gen from_sector = .
gen to_region = ""

* Generate Additional Names
foreach var in gexp dc dva vax davax ref ddc fc fva fdc gvc gvcb gvcf {
	gen `var' = .
}
	
foreach y in $years {
	
	di "Year: `y'"
	
	* Load ICIO Table
	icio_clean
	icio_load, iciot(user, userp($csv_path) tablen(${data}_`y'.csv) countrylist(${data}_countrylist.csv))  
	icio, info

	* Detailed Decomposition
	foreach ci in $countries {
		foreach ce in $countries {
			if "`ce'" == "`ci'" {
				continue	
			} 
			qui icio, exporter(`ce', all) importer(`ci') perspective(exporter) approach(source) output(detailed)
			// matlist r(detailed)
			mat res = r(detailed)
			local iter = `iter'+1
			local start = (`iter'-1)*$nsec + 1
			local end = `start' + $nsec - 1
			* di "EXP: `ce', IMP: `ci'-`si', Rows: `start'-`end'"
			local i = 0
			qui forvalues j = `start'/`end' {
				local i = `i'+1
				replace year = `y' in `j'
				replace from_region = "`ce'" in `j'
				replace from_sector = `i' in `j'
				replace to_region = "`ci'" in `j'
				replace gexp = res[1, `i'] in `j'
				replace dc = res[2, `i'] in `j'
				replace dva = res[3, `i'] in `j'
				replace vax = res[4, `i'] in `j'
				replace davax = res[5, `i'] in `j'
				replace ref = res[6, `i'] in `j'
				replace ddc = res[7, `i'] in `j'
				replace fc = res[8, `i'] in `j'
				replace fva = res[9, `i'] in `j'
				replace fdc = res[10, `i'] in `j'
				replace gvc = res[11, `i'] in `j'
				replace gvcb = res[12, `i'] in `j'
				replace gvcf = res[13, `i'] in `j'
			}			
		}
}
}


outsheet using "${csv_path}/${data}_GVC_BIL_SEC_BM19.csv", comma replace


