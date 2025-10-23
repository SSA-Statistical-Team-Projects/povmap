** uses something closer to the original EMDI code 
** point_estim has two slight modifications in nlme::lme, uses ML instead of REML fitting method and uses returnobject=true)
** initializes this in all parallel cores 
** relies on a specific path call where point_estimation.R, so will not run in other folders without modifying the path in mse_estimation.R

cap program drop Rpovmap
program define Rpovmap 
syntax namelist, smp_data(string) pop_data(string) smp_domains(string) pop_domains(string) [pop_subdomains(string) smp_subdomains(string) weights(string) WEIGHTS_Type(string) pop_weights(string) mse_pop_weights(string) threshold(string)  l(int 50) b(int 100) mse(string) transformation(string) interval(string) na_rm(string) cpus(int 1) seed(int 123) savexls(string) saveobject(string) interval(string) benchmark(string) aggregate_to(string) weights_type(string) benchmark_level(string) benchmark_type(string) benchmark_weights(string) rescale_weights(string) nlme_maxiter(int 1000) nlme_tolerance(real 1e-06) nlme_opt(string) boot_type(string) ydump(string) model_parameters(string) data_table(string) indicators(string) rpath(passthru)] 



* checks 
if real("`threshold'")==. & "`threshold'"~="" {
	noi dis "Threshold must be a real number or empty string"
	exit 
}
if "`saveobject'"=="" & "savexls"=="" {
	noi dis "Please specify either the savexls or saveobject option, and preferably both, to save your results"
	exit 
}
if "`benchmark'"~="" & "`benchmark_level'"=="" {
	noi dis "Please specify bencmhark_level option if benchmarking estimates"
	exit 
}


* set defaults 

if "`weights'"=="" {
local weights="NULL"
}
else {
	local weights `""`weights'""'
}

if "`weights_type'"=="" {
	local weights_type = "NULL"
}
else {
	local weights_type `""`weights_type'""'
}
if "`pop_weights'"=="" {
local pop_weights="NULL"
}
else {
	local pop_weights `""`pop_weights'""'
}

if "mse_pop_weights'"=="" {
	local mse_pop_weights="NULL" 
}
else {
	local mse_pop_weights `""`mse_pop_weights'""'
}

if "`threshold'"=="" {
local threshold="NULL"
}


if "`mse'"=="" {
local mse="FALSE"
}
if "`transformation'"=="" {
local transformation "box.cox"
}
if "`interval'"=="" {
local interval "default"
}
if "na_rm"=="" {
	local na_rm "FALSE"
} 
if "`savexls'"=="" {
local saving "EMDI_results.xlsx"
}
if "`saveobject'"=="" {
local saveobject "EMDI"
}
if "`scale_bootvar'"=="" {
local scale_bootvar "NULL"
}
if "`benchmark'"=="" {
	local benchmark "NULL"
}
else if strpos("`benchmark'","Head_Count")>0 | strpos("`benchmark'","Mean")>0 { 
	local benchmark `"c("`benchmark'")"'
}
else {
    local benchmark "`benchmark'"
	local benchmark_external "TRUE" 
} 

if "`benchmark_level'"=="" {
	local benchmark_level  "NULL"
}
else {
	local benchmark_level `""`benchmark_level'""'
}
if "`benchmark_type'"=="" {
	local benchmark_type  `""ratio""'  
}
else {
	local benchmark_type `""`benchmark_type'""'
}
if "`aggregate_to'"=="" {
	local aggregate_to "NULL"
}
else {
	local aggregate_to `""`aggregate_to'""'
}


if "`benchmark_weights'"=="" {
	local benchmark_weights  "NULL" 
}
else {
	local benchmark_weights `""`benchmark_weights'""'
}

if "`rescale_weights'"=="" {
	local rescale_weights  "FALSE"  
}
else {
	local rescale_weights "TRUE"
}
 
if "`boot_type'"=="" {
local boot_type="parametric"
}

if "`nlme_opt'"=="" {
local nlme_opt "nlminb"
}
 
  if "`pop_subdomains'"=="" {
  	local pop_subdomains "NULL"
  }
  else {
  	local pop_subdomains `""`pop_subdomains'""'
  }
 
 if "`smp_subdomains'"=="" {
  	local smp_subdomains "NULL"
  }
 else {
  	local smp_subdomains `""`smp_subdomains'""'
  }
  if `mse_pop_weights'=="" {
  	local mse_pop_weights "NULL"
  }
 else {
  	local mse_pop_weights `"`mse_pop_weights'"'
  }
  
if "`ydump'"!="" {
	local ydump `""`ydump'""'
}
else {
	local ydump "NULL"
} 


if "`model_parameters'"=="" {
	local model_parameters `""fixed""' 
}
else {
	local model_parameters `""model_parameters""' 
}

if "`data_table'"=="" {
	local data_table "FALSE"  
}

if "`indicators'"=="" {
	local vectorize "NULL"  
}
else {
	local indicators `"`indicators'"'
}





 
local pop_data : subinstr local pop_data "\" "/", all 
local smp_data : subinstr local smp_data "\" "/", all 
local savexls : subinstr local savexls "\" "/", all 
local saveobject : subinstr local saveobject "\" "/", all 

preserve
use `smp_data', replace 
* output model to `modelfile', a textfile R can read 
tempname modelfile  
gettoken yvar xvars : namelist     
local formula "`yvar' ~ "
foreach var of varlist `xvars' {
local formula "`formula' + `var'"
}
local formula : subinstr local formula "+ " ""
dis "`formula'"

cap file close R 
file open R using "`modelfile'.txt", write replace 
file write R "`formula'" _n
file close R
local working_dir : pwd
local working_dir : subinstr local working_dir "\" "/", all
restore 


file open Rscript using "povmap.R", write replace 
file write Rscript `"library("haven");"' _n
file write Rscript `"library("povmap");"' _n
file write Rscript `"setwd("`working_dir'");"' _n
file write Rscript `"pop <- as_factor(as.data.frame(read_dta("`pop_data'")));"' _n
file write Rscript `"smp <- as_factor(as.data.frame(read_dta("`smp_data'")));"' _n
file write Rscript `"model <- read.delim("`modelfile'.txt", header = FALSE, sep = "\t");"' _n
file write Rscript "model <- as.formula(as.character(model[1,1]));" _n

if "`benchmark_external'"=="TRUE" {
local benchmark_level_clean : list clean benchmark_level 	
file write Rscript "benchmark_value <- tapply(smp\$`benchmark', smp\$`benchmark_level_clean', mean)" _n     
file write Rscript "benchmark_table <- data.frame(`benchmark_level_clean' = names(benchmark_value), Head_Count = benchmark_value)" _n     
local benchmark "benchmark_table"
}


file write Rscript "ebp_results <- ebp(fixed = model,pop_data = pop," _n 
file write Rscript `"pop_domains = "`pop_domains'", smp_data = smp, smp_domains = "`smp_domains'","' _n 
file write Rscript `"pop_subdomains = `pop_subdomains', smp_subdomains = `smp_subdomains',"' _n 
file write Rscript `"threshold = `threshold', L = `l', B = `b', MSE = `mse', transformation = "`transformation'", interval = "`interval'","' _n 
file write Rscript `"boot_type="`boot_type'",na.rm = `na_rm', cpus = `cpus', seed=`seed', weights = `weights', weights_type = `weights_type', pop_weights = `pop_weights', MSE_pop_weights = `mse_pop_weights', aggregate_to = `aggregate_to',benchmark = `benchmark',"' _n          
file write Rscript `" benchmark_type = `benchmark_type', benchmark_level = `benchmark_level', benchmark_weights = `benchmark_weights', rescale_weights = `rescale_weights', nlme_maxiter = `nlme_maxiter',nlme_tolerance = `nlme_tolerance',nlme_opt = "`nlme_opt'",Ydump=`ydump',model_parameters=`model_parameters',data.table=`data_table',indicators=`indicators')"' _n          
if "`saveobject'"~="" {
file write Rscript `"save(ebp_results,file="`saveobject'")"' _n  	
}
file write Rscript `"write.excel(ebp_results, file = "`savexls'", indicator = "all", MSE = `mse', CV = `mse', split = FALSE,model=TRUE)"' _n          
file close Rscript 

rscript using povmap.R, require(povmap haven) `rpath'

end 
 
