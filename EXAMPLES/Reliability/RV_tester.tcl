# Random variable tester
# by Kevin Mackie 2012/04/12

set RV1 [list normal lognormal gamma shiftedExponential shiftedRayleigh exponential \
			  uniform type1LargestValue type1SmallestValue type2LargestValue \
			  chiSquare gumbel weibull laplace]
puts "RV set by mean and standard deviation: $RV1"
set nrv1 [llength $RV1]
set mean 5.0
set stdv 1.0

set RV2 [list beta rayleigh type3SmallestValue pareto]
puts "RV set by parameter only: $RV2"
set nrv2 [llength $RV2]

# create list of x-points for query
set xmin 0.0
set xmax 10.0
set npts 250
set dx [expr ($xmax-$xmin)/$npts]
set xpt [list $xmin]
for { set kl 1 } { $kl <= $npts } { incr kl } {
	lappend xpt [expr $xmin+$dx*$kl]
}

# create list of p-points for query
set pmin 1.0e-4
set pmax [expr 1.0-$pmin]
set dp [expr ($pmax-$pmin)/$npts]
set ppt [list $pmin]
for { set kl 1 } { $kl <= $npts } { incr kl } {
	lappend ppt [expr $pmin+$dp*$kl]
}

# CREATE THE RELIABILITY MODEL BUILDER ------------------------------------------------------------
reliability
model basic -ndm 2

# CREATE RANDOM VARIABLES
for { set kl 0 } { $kl < $nrv1 } { incr kl } {
	randomVariable  [expr $kl+1]	[lindex $RV1 $kl]		-mean $mean		-stdv $stdv
	parameter [expr $kl+1] randomVariable [expr $kl+1]
}

# TEST RVs
for { set kl 0 } { $kl < $nrv1 } { incr kl } {
	for { set jl 0 } { $jl <= $npts } { incr jl } {	
		set PDF($kl,$jl) [getPDF [expr $kl+1] [lindex $xpt $jl]]
		set CDF($kl,$jl) [getCDF [expr $kl+1] [lindex $xpt $jl]]
		set iCDF($kl,$jl) [getInverseCDF [expr $kl+1] [lindex $ppt $jl]]
	}
}

# dump to files for comparison with Matlab
for { set kl 0 } { $kl < $nrv1 } { incr kl } {
	set file [open "RV_tester_output/[lindex $RV1 $kl].out" "w"]
	for { set jl 0 } { $jl <= $npts } { incr jl } {	
		puts $file "[lindex $xpt $jl] $PDF($kl,$jl) $CDF($kl,$jl) [lindex $ppt $jl] $iCDF($kl,$jl)"
	}
	close $file
}

# do special RV individually that need parameters
wipeReliability
wipe
unset PDF CDF iCDF

reliability
model basic -ndm 2

randomVariable  1	[lindex $RV2 0]		-parameters 1 9 2 3
randomVariable  2	[lindex $RV2 1]		-parameters [expr 2.0*($mean*2.0/3.0)/sqrt(acos(-1.0))]
randomVariable  3	[lindex $RV2 2]		-parameters 0.5 3 4
randomVariable  4	[lindex $RV2 3]		-parameters 3 4

for { set kl 0 } { $kl < $nrv2 } { incr kl } {
	parameter [expr $kl+1] randomVariable [expr $kl+1]
}

# TEST RVs
for { set kl 0 } { $kl < $nrv2 } { incr kl } {
	for { set jl 0 } { $jl <= $npts } { incr jl } {	
		set PDF($kl,$jl) [getPDF [expr $kl+1] [lindex $xpt $jl]]
		set CDF($kl,$jl) [getCDF [expr $kl+1] [lindex $xpt $jl]]
		set iCDF($kl,$jl) [getInverseCDF [expr $kl+1] [lindex $ppt $jl]]
	}
}

# dump to files for comparison with Matlab
for { set kl 0 } { $kl < $nrv2 } { incr kl } {
	set file [open "RV_tester_output/[lindex $RV2 $kl].out" "w"]
	for { set jl 0 } { $jl <= $npts } { incr jl } {	
		puts $file "[lindex $xpt $jl] $PDF($kl,$jl) $CDF($kl,$jl) [lindex $ppt $jl] $iCDF($kl,$jl)"
	}
	close $file
}

# test incomplete gamma function using the gamma CDF
wipeReliability
wipe
unset CDF PDF iCDF

reliability
model basic -ndm 2

set klist [list 0.1 0.25 0.5 0.75 1.0 1.25 2.5 5.0 10.0 25.0]
set nks [llength $klist]

# create the variables
for { set kl 0 } { $kl < $nks } { incr kl } {
	randomVariable  [expr $kl+1]	gamma		-parameters [lindex $klist $kl] 1.0
	parameter [expr $kl+1] randomVariable [expr $kl+1]
}

# TEST RVs
for { set kl 0 } { $kl < $nks } { incr kl } {
	for { set jl 0 } { $jl <= $npts } { incr jl } {	
		set PDF($kl,$jl) [getPDF [expr $kl+1] [lindex $xpt $jl]]
		set CDF($kl,$jl) [getCDF [expr $kl+1] [lindex $xpt $jl]]
	}
}

# dump to files for comparison with Matlab
for { set kl 0 } { $kl < $nks } { incr kl } {
	set file [open "RV_tester_output/gamma[expr $kl+1].out" "w"]
	for { set jl 0 } { $jl <= $npts } { incr jl } {	
		puts $file "[lindex $xpt $jl] $PDF($kl,$jl) $CDF($kl,$jl)"
	}
	close $file
}
