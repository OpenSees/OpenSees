###########################################################
#                                                         #
# Procedure to compute ultimate resistance, tult, and     #
#  displacement at 50% mobilization of tult, z50, for     #
#  use in t-z curves for cohesionless soil.               #
#                                                         #
#   Created by:  Chris McGann                             #
#                University of Washington                 #
#                                                         #
###########################################################

proc get_tzParam { phi b sigV pEleLength} {

# references
#  Mosher, R.L. (1984). “Load transfer criteria for numerical analysis of
#   axial loaded piles in sand.” U.S. Army Engineering and Waterways
#   Experimental Station, Automatic Data Processing Center, Vicksburg, Miss.
#
#  Kulhawy, F.H. (1991). "Drilled shaft foundations." Foundation engineering
#   handbook, 2nd Ed., Chap 14, H.-Y. Fang ed., Van Nostrand Reinhold, New York

    set pi 3.14159265358979
    
  # Compute tult based on tult = Ko*sigV*pi*dia*tan(delta), where
  #   Ko    is coeff. of lateral earth pressure at rest, 
  #         taken as Ko = 0.4
  #   delta is interface friction between soil and pile,
  #         taken as delta = 0.8*phi to be representative of a 
  #         smooth precast concrete pile after Kulhawy (1991)
    set delta [expr 0.8*$phi*$pi/180]

  # if z = 0 (ground surface) need to specify a small non-zero value of sigV
    if { $sigV == 0.0 } {
        set sigV 0.01
    }

    set tu   [expr 0.4*$sigV*$pi*$b*tan($delta)]
    
  # TzSimple1 material formulated with tult as force, not stress, multiply by tributary length of pile
    set tult [expr $tu*$pEleLength]

  # Mosher (1984) provides recommended initial tangents based on friction angle
	# values are in units of psf/in
    set kf(1)    6000 
    set kf(2)    10000
	set kf(3)    10000
	set kf(4)    14000
	set kf(5)    14000
	set kf(6)    18000

    set fric(1)  28
    set fric(2)  31
    set fric(3)  32
    set fric(4)  34
    set fric(5)  35
    set fric(6)  38

    set dataNum  6

	# determine kf for input value of phi, linear interpolation for intermediate values
    if { $phi < $fric(1) } {
        set k $kf(1)
    } elseif { $phi > $fric(6) } {
        set k $kf(6)
    } else {
        for {set i 1} {$i <= [expr $dataNum-1] } {incr i} {
            if { ($fric($i) <= $phi) && ($phi <= $fric([expr $i+1])) } {
                set k [expr ($kf([expr $i+1]) - $kf($i))/($fric([expr $i+1]) - $fric($i))*($phi - $fric($i)) + $kf($i)]
            } 
        }
    }

  # need to convert kf to units of kN/m^3
    set kSIunits [expr $k*1.885]

  # based on a t-z curve of the shape recommended by Mosher (1984), z50 = tult/kf
    set z50 [expr $tult/$kSIunits]

  # return values of tult and z50 for use in t-z material
    set outResult [concat $tult $z50]

    return $outResult
}
