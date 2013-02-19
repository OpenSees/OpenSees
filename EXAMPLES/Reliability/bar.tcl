model basic -ndm 1

reliability

set R 2000.0
set S 1700.0

randomVariable 1 normal -mean $R -stdv 135.0
randomVariable 2 normal -mean $S -stdv 200.0

parameter 1 randomVariable 1
parameter 2 randomVariable 2

performanceFunction 1 "\$par(1)-\$par(2)"
performanceFunction 2 "log(\$par(1))-log(\$par(2))"
performanceFunction 3 "pow(\$par(1),2)-pow(\$par(2),2)"

functionEvaluator Tcl
gradientEvaluator FiniteDifference -pert 1000
startPoint Mean

runFOSMAnalysis barFOSM.out
