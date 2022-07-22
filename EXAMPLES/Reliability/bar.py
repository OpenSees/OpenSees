import opensees as ops

ops.wipe()

ops.model('basic', '-ndm', 1)


R = 2000.0

S = 1700.0


ops.randomVariable(1, 'normal', '-mean', R, '-stdv', 135.0)

ops.randomVariable(2, 'normal', '-mean', S, '-stdv', 200.0)


ops.parameter(1, 'randomVariable', 1)

ops.parameter(2, 'randomVariable', 2)


ops.performanceFunction(1, "par[1]-par[2]")

ops.performanceFunction(2, "log(par[1])-log(par[2])")

ops.performanceFunction(3, "pow(par[1],2)-pow(par[2],2)")

ops.functionEvaluator('Python')

ops.gradientEvaluator('FiniteDifference', '-pert', 1000)

ops.startPoint('Mean')

ops.randomNumberGenerator('CStdLib')

ops.reliabilityConvergenceCheck('Standard')

ops.searchDirection('iHLRF')

ops.meritFunctionCheck('AdkZhang')

ops.probabilityTransformation('Nataf')

ops.stepSizeRule('Armijo')

ops.rootFinding('Secant')



ops.runFOSMAnalysis('barFOSM.out')
