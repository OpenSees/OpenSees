# source for opstkUniaxialMaterialViewer.tcl
# define sample input, add to Material menu

$m add command -label "Elastic" -command "uniaxialMaterialSample Elastic matTag=1 E=25."
$m add command -label "ElasticPP" -command "uniaxialMaterialSample ElasticPP matTag=1 E=100. yieldStrain=0.25"
$m add command -label "Steel01" -command "uniaxialMaterialSample Steel01 matTag=1 Fy=40. E=100. b=0.1"
$m add command -label "Steel02" -command "uniaxialMaterialSample Steel02 matTag=1 Fy=40. E=100. b=0.1 R0=15 CR1=0.925 CR2=0.15"
$m add command -label "MultiLinear" -command "uniaxialMaterialSample MultiLinear matTag=1 e1=0.05 s1=10 e2=0.15 s2=20 e3=0.25 s3=25 e4=0.6 s4=0."
$m add command -label "SimpleFracture" -command "uniaxialMaterialSample SimpleFracture matTag=1 otherTag=1 maxStrain=0.4"
