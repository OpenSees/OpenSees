import opensees as ops

ops.model('basic', '-ndm', 2, '-ndf', 2)

L = 144.0

ops.node(1, 0, 0)
ops.node(2, L, 0.0)

ops.fix(1, 1, 1)
ops.fix(2, 0, 1)

E = 30000.0
A = 25.0
fy = 50.0

ops.uniaxialMaterial("Hardening", 1, E, fy, 0, 100.0)

ops.element("truss", 1, 1, 2, A, 1)

P = 25.0

tsTag = 1
ops.timeSeries("Linear", tsTag)

patternTag = 1
ops.pattern("Plain", patternTag, tsTag)

ops.load(2, P, 0)

ops.analysis("Static")

ops.randomVariable(62, 'lognormal', '-mean', E,
                   '-stdv', 0.1*E)
ops.randomVariable(25, 'lognormal', '-mean', A,
                   '-stdv', 0.1*A)
ops.randomVariable(33, 'lognormal', '-mean', fy,
                   '-stdv', 0.1*fy)

ops.parameter(12, 'randomVariable', 62,
              "element", 1, "E")
ops.parameter(13, 'randomVariable', 25,
              "element", 1, "A")
ops.parameter(14, 'randomVariable', 33,
              "element", 1, "fy")
ops.parameter(23, "node", 2, "disp", 1)

ops.performanceFunction(76, "5.5-par[23]")

ops.sensitivityAlgorithm("-computeAtEachStep")

ops.randomNumberGenerator('CStdLib')
ops.probabilityTransformation('Nataf', "-print", 1)
ops.reliabilityConvergenceCheck('Standard', "-e1", 1e-2,
                                "-e2", 1e-2, "-print", 1)
ops.functionEvaluator('Python', "-file", "opensees.analyze(55)")
ops.gradientEvaluator('Implicit')
ops.searchDirection('iHLRF')
ops.meritFunctionCheck('AdkZhang', "-multi", 2.0, "-add", 10.0,
                       "-factor", 0.5)
ops.stepSizeRule('Fixed', "-stepSize", 1.0)
ops.startPoint('Mean')
ops.findDesignPoint('StepSearch', "-maxNumIter", 15,
                    "-printDesignPointX", "designPointX.out")

ops.runFORMAnalysis('py_truss_FORM.out')

for perf in ops.getLSFTags():
  print(f"Performance Function {perf}")
  print(f"beta = {ops.betaFORM[perf]:.4f}")

ops.runFOSMAnalysis('py_truss_FOSM.out')

# Now reset probability transformation so we don't get as much output for simulation
ops.probabilityTransformation('Nataf', "-print", 0)

ops.runImportanceSamplingAnalysis("py_truss_SIM1.out",
                                  "-type", "failureProbability",
                                  "-variance", 1.0, "-maxNum", 5000,
                                  "-targetCOV", 0.01, "-print", 0)

ops.runImportanceSamplingAnalysis("py_truss_SIM2.out",
                                  "-type", "responseStatistics",
                                  "-variance", 1.0, "-maxNum", 5000,
                                  "-targetCOV", 0.01, "-print", 0)
