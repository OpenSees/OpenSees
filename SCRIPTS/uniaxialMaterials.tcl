# source for opstkUniaxialMaterialViewer.tcl
# define sample input, add to Material menu

$m add command -label "Elastic" -command "uniaxialMaterialSample Elastic matTag=1 E=25."
$m add command -label "ElasticPP" -command "uniaxialMaterialSample ElasticPP matTag=1 E=100. yieldStrain=0.25"
$m add command -label "Steel01" -command "uniaxialMaterialSample Steel01 matTag=1 Fy=40. E=100. b=0.1"
$m add command -label "Steel02" -command "uniaxialMaterialSample Steel02 matTag=1 Fy=40. E=100. b=0.1 R0=15 CR1=0.925 CR2=0.15"
$m add command -label "MultiLinear" -command "uniaxialMaterialSample MultiLinear matTag=1 e1=0.05 s1=10 e2=0.15 s2=20 e3=0.25 s3=25 e4=0.6 s4=0."
$m add command -label "SimpleFracture" -command "uniaxialMaterialSample SimpleFracture matTag=1 otherTag=1 maxStrain=0.4"
$m add command -label "MultiLinear2" -command "uniaxialMaterialSample MultiLinear matTag=1 e1=0.0001341737342935323 s1=24640.0 e2=0.008252161648910277 s2=25872.000000000004 e3=0.022357133345899175 s3=0.0 e4=10.0223571333459 s4=0."
$m add command -label "Bilin02" -command "uniaxialMaterialSample Bilin02 matTag=1 K0=3672850.0 asP=0.00082 asN=0.00082 myP=24640.0 myN=-24640.0 LS=1.60688 LC=1.60688 LA=1.60688 LK=1.60688 cS=1.0 cC=1.0 cA=1.0 cK=1.0 tpP=0.405899 tpN=0.405899 tpcP=0.705248 tpcN=0.705248 rP=0.4 rN=0.4 tuP=0.4 tuN=0.4 dP=1.0 dN=1.0 nFactor=1.";
$m add command -label "ModIMKPeak" -command "uniaxialMaterialSample ModIMKPeakOriented matTag=1 K0=3672850.0 asP=0.00082 asN=0.00082 myP=24640.0 myN=-24640.0 LS=1.60688 LC=1.60688 LA=1.60688 LK=1.60688 cS=1.0 cC=1.0 cA=1.0 cK=1.0 tpP=0.405899 tpN=0.405899 tpcP=0.705248 tpcN=0.705248 rP=0.4 rN=0.4 tuP=0.4 tuN=0.4 dP=1.0 dN=1.0";
$m add command -label "Bilin" -command "uniaxialMaterialSample Bilin matTag=1 p1=3672850.0 p2=0.0008263977213611486 p3=0.0008263977213611486 p4=24640.0 p5=-24640.0 p6=1.60688 p7=1.60688 p8=1.6068826060003187 p9=1.60688 p10=1.0 p11=1.0 p12=1.0 p13=1.0 p14=0.405899 p15=0.405899 p16=0.705248 p17=0.705248 p18=0.4 p19=0.4 p20=0.4 p21=0.4 p22=1.0 p23=1.0 p24=1.";
$m add command -label "HoehlerStanton" -command "uniaxialMaterialSample HoehlerStanton matTag=1 200000 0.015 400 550 0.15 3 25 23 0.2"
$m add command -label "BoucWenOriginal" -command "uniaxialMaterialSample BoucWenOriginal matTag=1 E=41760.0 fy=40.0 alphaL=0.025 alphaNL=0.0 mu=2.0 eta=1.0 beta=0.55 gamma=0.45"
$m add command -label "BoucWen" -command "uniaxialMaterialSample BoucWen matTag=1 alpha=0.025 ko=41760.0 n=1.0 gamma=0.45 beta=0.55 Ao=1.0 deltaA=0.0 deltaNu=0.0 deltaEta=0.0"