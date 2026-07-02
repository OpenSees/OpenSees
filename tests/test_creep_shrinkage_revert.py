try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops
import random
from statistics import NormalDist
from math import isclose

_std_normal = NormalDist()


def _modulator(t, t1, t2, t3, t4):
    if t <= t1 or t >= t4:
        return 0.0
    elif t < t2:
        return (t - t1) / (t2 - t1)
    elif t < t3:
        return 1.0
    else:
        return 1.0 - (t - t3) / (t4 - t3)


def white_noise(dt, Tfinal, tup, tdown, seed=None):
    """Trapezoidal-modulated white noise strain history.
    Adapted from bennycloth.timeseries.whitenoise.WhiteNoise (Michael H.
    Scott) -- inlined here to avoid bennycloth's package-level import of
    openseespy, which is unrelated to WhiteNoise itself and pulls in a
    broken PyPI Linux binary as an unnecessary dependency."""
    if seed is not None:
        random.seed(seed)

    Npts = int(Tfinal / dt)
    noise = [_std_normal.inv_cdf(random.random()) for _ in range(Npts)]

    t1, t2, t3, t4 = 0, tup, Tfinal - tdown, Tfinal
    modnoise = [noise[i] * _modulator(dt * (i + 1), t1, t2, t3, t4) for i in range(Npts)]

    maxnoise = max(abs(v) for v in modnoise)
    return [v / maxnoise for v in modnoise]

# Base material is Concrete02IS, matching the CreepShrinkageACI209
# documentation example (7-day moist-cured concrete loaded at 28 days).
# CreepShrinkageACI209 is a concrete creep/shrinkage wrapper, so a concrete
# base material is required to exercise realistic (compression-dominated,
# nonlinear) behavior -- a steel material would not exercise the intended
# use case.
E0 = 4000.0
fpc = -4.0
epsc0 = -0.002
fpcu = -0.8
epscu = -0.01

dtp = 1.0    # time between white-noise pulses
T = 40.0     # total time series duration
t_up = T / 4
t_down = T / 10
eps_max = 0.5 * abs(epsc0)  # keep well inside the linear/pre-peak range

dt = dtp / 20.0  # analysis time step
Nsteps = int(T / dt)


def define_model(creep_on):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0); ops.fix(1, 1)
    ops.node(2, 0)

    ops.setCreep(1 if creep_on else 0)

    ops.uniaxialMaterial('Concrete02IS', 1, E0, fpc, epsc0, fpcu, epscu)
    ops.uniaxialMaterial(
        'CreepShrinkageACI209', 2, 1,
        7.0,      # tD
        -780e-6,  # epsshu
        35.0,     # psish
        28.0,     # Tcr
        2.35,     # phiu
        0.6,      # psicr1
        10.0,     # psicr2
        0.0,      # tcast
    )

    ops.element('zeroLength', 1, 1, 2, '-mat', 2, '-dir', 1)

    strain_history = white_noise(dtp, T, t_up, t_down, seed=221537)

    ops.timeSeries('Path', 1, '-dt', dtp, '-values', *strain_history, '-factor', eps_max)
    ops.pattern('Plain', 1, 1)
    ops.sp(2, 1, 1.0)

    ops.integrator('LoadControl', dt)
    ops.constraints('Transformation')
    ops.system('UmfPack')
    ops.algorithm('Newton')
    ops.analysis('Static', '-noWarnings')


def get_baseline(creep_on):
    define_model(creep_on)
    ops.test('NormUnbalance', 1e-8, 10, 0)

    sig0 = [0.0] * Nsteps
    for i in range(Nsteps):
        ok = ops.analyze(1)
        assert ok == 0, f"baseline analysis failed to converge at step {i}"
        sig0[i] = ops.eleResponse(1, 'material', 1, 'stress')[0]
    return sig0


def check_revert_to_last_commit(creep_on):
    """Force every step to fail once (negative tolerance -> revertToLastCommit),
    then retry with a positive tolerance. Stress history must match baseline."""
    sig0 = get_baseline(creep_on)

    define_model(creep_on)
    for i in range(Nsteps):
        ops.test('NormUnbalance', -1e-8, 10, 0)
        ops.analyze(1)  # expected to fail; triggers revertToLastCommit

        ops.test('NormUnbalance', 1e-8, 10, 0)
        ok = ops.analyze(1)
        assert ok == 0, f"retry analysis failed to converge at step {i}"

        sigi = ops.eleResponse(1, 'material', 1, 'stress')[0]
        assert isclose(sigi, sig0[i], abs_tol=1e-9, rel_tol=1e-9), \
            f"revertToLastCommit diverged at step {i} (t={ (i+1)*dt }): {sigi} vs {sig0[i]}"


def check_revert_to_start(creep_on, repeats=3):
    """Repeat the baseline analysis (no forced failures) after ops.reset(),
    without ops.wipe(). All repeats must produce identical results."""
    define_model(creep_on)
    ops.test('NormUnbalance', 1e-8, 10, 0)

    reference = None
    for r in range(repeats):
        ops.reset()
        sig = [0.0] * Nsteps
        for i in range(Nsteps):
            ok = ops.analyze(1)
            assert ok == 0, f"repeat {r} step {i} failed to converge"
            sig[i] = ops.eleResponse(1, 'material', 1, 'stress')[0]

        if reference is None:
            reference = sig
        else:
            for i in range(Nsteps):
                assert isclose(sig[i], reference[i], abs_tol=1e-9, rel_tol=1e-9), \
                    f"revertToStart (reset) diverged on repeat {r}, step {i}: {sig[i]} vs {reference[i]}"


def test_revert_to_last_commit_creep_on():
    check_revert_to_last_commit(True)


def test_revert_to_last_commit_creep_off():
    check_revert_to_last_commit(False)


def test_revert_to_start_creep_on():
    check_revert_to_start(True)


def check_revert_multi_fail(creep_on, num_fails=3):
    """Force several consecutive failed attempts per step (each one triggering
    its own revertToLastCommit) before the final successful retry."""
    sig0 = get_baseline(creep_on)

    define_model(creep_on)
    for i in range(Nsteps):
        for _ in range(num_fails):
            ops.test('NormUnbalance', -1e-8, 10, 0)
            ops.analyze(1)  # expected to fail; triggers revertToLastCommit

        ops.test('NormUnbalance', 1e-8, 10, 0)
        ok = ops.analyze(1)
        assert ok == 0, f"retry analysis failed to converge at step {i}"

        sigi = ops.eleResponse(1, 'material', 1, 'stress')[0]
        assert isclose(sigi, sig0[i], abs_tol=1e-9, rel_tol=1e-9), \
            f"revertToLastCommit (multi-fail) diverged at step {i} (t={ (i+1)*dt }): {sigi} vs {sig0[i]}"


def check_revert_with_algorithm_switch(creep_on):
    """Mimic 'smart analyze' logic: fail with Newton, switch to a different
    algorithm to retry (and switch back), matching the exact scenario the
    blog post names as the real-world trigger for revertToLastCommit()."""
    sig0 = get_baseline(creep_on)

    define_model(creep_on)
    for i in range(Nsteps):
        ops.algorithm('Newton')
        ops.test('NormUnbalance', -1e-8, 10, 0)
        ops.analyze(1)  # expected to fail; triggers revertToLastCommit

        ops.algorithm('ModifiedNewton')
        ops.test('NormUnbalance', 1e-8, 10, 0)
        ok = ops.analyze(1)
        assert ok == 0, f"retry (ModifiedNewton) failed to converge at step {i}"

        ops.algorithm('Newton')  # restore for next step's baseline-matching attempt

        sigi = ops.eleResponse(1, 'material', 1, 'stress')[0]
        assert isclose(sigi, sig0[i], abs_tol=1e-9, rel_tol=1e-9), \
            f"revertToLastCommit (algorithm switch) diverged at step {i} (t={ (i+1)*dt }): {sigi} vs {sig0[i]}"


def test_revert_to_start_creep_off():
    check_revert_to_start(False)


def test_revert_multi_fail_creep_on():
    check_revert_multi_fail(True)


def test_revert_multi_fail_creep_off():
    check_revert_multi_fail(False)


def test_revert_with_algorithm_switch_creep_on():
    check_revert_with_algorithm_switch(True)


def test_revert_with_algorithm_switch_creep_off():
    check_revert_with_algorithm_switch(False)


# ---------------------------------------------------------------------------
# Reproduction attempt matching openseesdigital.com/white-noise/uniaxial/
# Material 102, which reports:
#   - PASS reset()
#   - FAIL revertToLastCommit() at timestep 651 (max diff ~0.0044)
#   - FAIL restore() at timestep 21 (max diff ~0.0044)
# Parameters below are taken directly from that page's material definition:
#   ops.uniaxialMaterial('Concrete01', 23, -30, -0.002, -6, 0.006)
#   ops.uniaxialMaterial('CreepShrinkageACI209', 1, 23, 2, -600e-6, 42.0,
#                         100, 3.0, 1, 42.0, 0)
# i.e. a much larger phiu/Tcr/psish than the docs example, over many more
# steps and a wider strain range that drives well past epsc0/epscu into the
# post-peak/softening branch and into tension.
# ---------------------------------------------------------------------------

W_dtp = 1.0
W_T = 250.0
W_t_up = W_T / 4
W_t_down = W_T / 10
W_dt = W_dtp / 20.0
W_Nsteps = int(W_T / W_dt)
W_eps_max = 0.006  # matches the ~0.006-0.015 strain range implied by the page's plots


def define_model_website(creep_on):
    ops.wipe()
    ops.model('basic', '-ndm', 1, '-ndf', 1)
    ops.node(1, 0); ops.fix(1, 1)
    ops.node(2, 0)

    ops.setCreep(1 if creep_on else 0)

    ops.uniaxialMaterial('Concrete01', 23, -30, -0.002, -6, 0.006)
    ops.uniaxialMaterial('CreepShrinkageACI209', 1, 23, 2, -600e-6, 42.0, 100, 3.0, 1, 42.0, 0)

    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)

    strain_history = white_noise(W_dtp, W_T, W_t_up, W_t_down, seed=221537)

    ops.timeSeries('Path', 1, '-dt', W_dtp, '-values', *strain_history, '-factor', W_eps_max)
    ops.pattern('Plain', 1, 1)
    ops.sp(2, 1, 1.0)

    ops.integrator('LoadControl', W_dt)
    ops.constraints('Transformation')
    ops.system('UmfPack')
    ops.algorithm('Newton')
    ops.analysis('Static', '-noWarnings')


def get_baseline_website(creep_on):
    define_model_website(creep_on)
    ops.test('NormUnbalance', 1e-8, 10, 0)

    sig0 = [0.0] * W_Nsteps
    for i in range(W_Nsteps):
        ok = ops.analyze(1)
        assert ok == 0, f"baseline analysis failed to converge at step {i}"
        sig0[i] = ops.eleResponse(1, 'material', 1, 'stress')[0]
    return sig0


def test_revert_to_last_commit_website_params_creep_on():
    sig0 = get_baseline_website(True)

    define_model_website(True)
    for i in range(W_Nsteps):
        ops.test('NormUnbalance', -1e-8, 10, 0)
        ops.analyze(1)  # expected to fail; triggers revertToLastCommit

        ops.test('NormUnbalance', 1e-8, 10, 0)
        ok = ops.analyze(1)
        assert ok == 0, f"retry analysis failed to converge at step {i}"

        sigi = ops.eleResponse(1, 'material', 1, 'stress')[0]
        assert isclose(sigi, sig0[i], abs_tol=1e-9, rel_tol=1e-9), \
            f"revertToLastCommit (website params) diverged at step {i} (t={ (i+1)*W_dt }): {sigi} vs {sig0[i]}"


def test_revert_to_start_website_params_creep_on(repeats=3):
    define_model_website(True)
    ops.test('NormUnbalance', 1e-8, 10, 0)

    reference = None
    for r in range(repeats):
        ops.reset()
        sig = [0.0] * W_Nsteps
        for i in range(W_Nsteps):
            ok = ops.analyze(1)
            assert ok == 0, f"repeat {r} step {i} failed to converge"
            sig[i] = ops.eleResponse(1, 'material', 1, 'stress')[0]

        if reference is None:
            reference = sig
        else:
            for i in range(W_Nsteps):
                assert isclose(sig[i], reference[i], abs_tol=1e-9, rel_tol=1e-9), \
                    f"revertToStart (website params) diverged on repeat {r}, step {i}: {sig[i]} vs {reference[i]}"
