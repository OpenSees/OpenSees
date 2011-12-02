#================================================================
# CreateMaterialModels.tcl
# --Create the all elements from file elements.dat
# Zhaohui Yang UC Davis
# March 31, 2002
#================================================================

# Elastic material model
set Eb  [expr 17400*$kPa]; #Young's Modulus or any number
set Vb 0.3; #Poisson's Ratio
set rhob [expr 1650.0*$kg/pow($m, 3)]; #unit weight

set Eo  [expr 14000*$kPa]; 
set Vo 0.3;  
set rhoo [expr 1500.0*$kg/pow($m, 3)]; 

set Ei  [expr 15000*$kPa]; 
set Vi 0.3;  
set rhoi [expr 1500.0*$kg/pow($m, 3)]; 

#================================================================
# Create an elastic material for testing
nDMaterial ElasticIsotropic3D  1   $Ei  $Vi  $rhoi
nDMaterial ElasticIsotropic3D  2   $Eb  $Vb  $rhob
nDMaterial ElasticIsotropic3D  3   $Eo  $Vo  $rhoo



