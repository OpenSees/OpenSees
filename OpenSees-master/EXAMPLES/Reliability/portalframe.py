import opensees as ops

ops.model('basic', '-ndm', 2, '-ndf', 3)

ops.node(1, 0, 0)
ops.node(2, 0, 144)
ops.node(3, 240, 144)
ops.node(4, 240, 0)

ops.fix(1, 1, 1, 1)
ops.fix(4, 1, 1, 1)

E = 30000.0
Ag = 25.0
Ig = 1500.0
Ac = 29.0
Ic = 2000.0

gsecTag = 1
ops.section("Elastic", gsecTag, E, Ag, Ig)

csecTag = 2
ops.section("Elastic", csecTag, E, Ac, Ic)

transfTag = 1
ops.geomTransf("Linear", transfTag)

N = 3

gbiTag = 1
ops.beamIntegration("Lobatto", gbiTag, gsecTag, N)
cbiTag = 2
ops.beamIntegration("Lobatto", cbiTag, csecTag, N)

leftColTag = 1
ops.element("forceBeamColumn", leftColTag, 1, 2, transfTag, cbiTag)
girderTag = 2
ops.element("forceBeamColumn", girderTag, 2, 3, transfTag, gbiTag)
rightColTag = 3
ops.element("forceBeamColumn", rightColTag, 3, 4, transfTag, cbiTag)

P = 25.0
w = 1.0e-1

tsTag = 1
ops.timeSeries("Constant", tsTag)

patternTag = 1
ops.pattern("Plain", patternTag, tsTag)

ops.load(2, P, 0, 0)
ops.eleLoad("-ele", girderTag, "-type", "beamUniform", -w)

ops.analysis("Static")

ops.randomVariable(62, 'lognormal', '-mean', E, '-stdv', 0.1*E)
ops.randomVariable(32, 'normal', '-mean', P, '-stdv', 0.2*P)
ops.randomVariable(89, 'normal', '-mean', 0, '-stdv', 1)
ops.randomVariable(41, 'normal', '-mean', -w, '-stdv', abs(0.2*w))

ops.parameter(12, 'randomVariable', 62,
              "element", leftColTag, "E")
ops.addToParameter(12, 'element', rightColTag, "E")
ops.parameter(25, 'randomVariable', 32,
              "loadPattern", patternTag,
              "loadAtNode", 2, 1)
ops.parameter(3, 'randomVariable', 89,
              "node", 1, "coord", 1)
ops.parameter(45, 'randomVariable', 41,
              "loadPattern", patternTag,
              "elementLoad", girderTag, "wy")
ops.parameter(23, "node", 2, "disp", 1)

ops.performanceFunction(76, "0.15-par[23]")

ops.sensitivityAlgorithm("-computeAtEachStep")

ops.randomNumberGenerator('CStdLib')
ops.probabilityTransformation('Nataf', "-print", 3)
ops.reliabilityConvergenceCheck('Standard', "-e1", 1e-2,
                                "-e2", 1e-2, "-print", 1)
ops.functionEvaluator('Python', "-file", "opensees.analyze(1)")
ops.gradientEvaluator('Implicit')
ops.searchDirection('iHLRF')
ops.meritFunctionCheck('AdkZhang', "-multi", 2.0, "-add", 50.0,
                       "-factor", 0.5)
ops.stepSizeRule('Armijo', "-maxNum", 10, "-base", 0.5,
                 "-initial", 0.3, 5)
ops.stepSizeRule('Fixed', "-stepSize", 1.0)
ops.startPoint('Mean')
ops.findDesignPoint('StepSearch', "-maxNumIter", 30)

ops.runFORMAnalysis('portalframe.out')

for perf in ops.getLSFTags():
    print(f"Performance Function {perf}")
    print(f"beta = {ops.betaFORM[perf]:.7f}")
    for rv in ops.getRVTags():
        print(
            f"\t x*({rv}) = {ops.designPointXFORM[perf, rv]:7.4f}, alpha({rv}) = {ops.alphaFORM[perf, rv]:7.4f}, gamma({rv}) = {ops.gammaFORM[perf, rv]}")
