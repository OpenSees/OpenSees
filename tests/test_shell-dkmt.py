"""
Tests for the ShellDKMT element (CST membrane + DKMT plate bending).

1. The 18 x 18 stiffness of one element in a general 3D position matches an
   independent pure-Python implementation of the formulation.
2. Rigid body motions are stress free (the drilling penalty is objective).
3. A thin simply supported plate gives the same deflection as ShellDKGT
   (DKMT reduces to DKT) and the Kirchhoff (Navier) solution.
4. A thick simply supported plate (a/h = 10) gives the Mindlin (FSDT)
   solution.
5. The self weight adds up to rho*h*A.
"""

import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# ---------------------------------------------------------------------------
# pure-Python reference of the formulation
# ---------------------------------------------------------------------------

KAPPA = 5.0 / 6.0


def _sub(a, b):
    return [a[i] - b[i] for i in range(len(a))]


def _dot(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _unit(a):
    n = math.sqrt(_dot(a, a))
    return [x / n for x in a]


def _zeros(n, m):
    return [[0.0] * m for _ in range(n)]


def _matmul(A, B):
    Bt = list(zip(*B))
    return [[_dot(row, col) for col in Bt] for row in A]


def _transpose(A):
    return [list(r) for r in zip(*A)]


def _solve(A, B):
    """Gauss elimination with partial pivoting, A n x n, B n x m"""
    n = len(A)
    M = [list(A[i]) + list(B[i]) for i in range(n)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(len(M[r]))]
    return [[M[i][n + k] / M[i][i] for k in range(len(B[0]))] for i in range(n)]


def ref_stiffness(X, E, nu, h, drill=1.0e-4):
    """global 18 x 18 stiffness, nodal dofs u v w rx ry rz"""
    e1 = _unit(_sub(X[1], X[0]))
    e3 = _unit(_cross(_sub(X[1], X[0]), _sub(X[2], X[0])))
    e2 = _cross(e3, e1)
    R = [e1, e2, e3]
    x = [[_dot(_sub(Xi, X[0]), e1), _dot(_sub(Xi, X[0]), e2)] for Xi in X]

    J = [_sub(x[1], x[0]), _sub(x[2], x[0])]
    det = J[0][0] * J[1][1] - J[0][1] * J[1][0]
    A = 0.5 * det
    Ji = [[J[1][1] / det, -J[0][1] / det], [-J[1][0] / det, J[0][0] / det]]
    dNxi = [[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]]
    dN = [[Ji[0][0] * d[0] + Ji[0][1] * d[1], Ji[1][0] * d[0] + Ji[1][1] * d[1]] for d in dNxi]

    G = E / (2.0 * (1.0 + nu))
    Dm = [[E * h / (1 - nu ** 2) * v for v in row] for row in [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]]]
    Db_ = E * h ** 3 / (12 * (1 - nu ** 2))
    Db = [[Db_ * v for v in row] for row in [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]]]
    Ds = KAPPA * G * h

    K = _zeros(18, 18)

    # membrane, dofs (u, v) of each node
    im = [6 * i + c for i in range(3) for c in (0, 1)]
    Bm = _zeros(3, 6)
    for i in range(3):
        Bm[0][2 * i] = dN[i][0]
        Bm[1][2 * i + 1] = dN[i][1]
        Bm[2][2 * i] = dN[i][1]
        Bm[2][2 * i + 1] = dN[i][0]
    Km = _matmul(_transpose(Bm), _matmul(Dm, Bm))
    for a in range(6):
        for b in range(6):
            K[im[a]][im[b]] += A * Km[a][b]

    # bending, dofs (w, beta_x, beta_y) of each node, beta_x = ry, beta_y = -rx
    edges = ((0, 1), (1, 2), (2, 0))
    C, S, L = [], [], []
    for i, j in edges:
        d = _sub(x[j], x[i])
        l = math.sqrt(_dot(d, d))
        C.append(d[0] / l)
        S.append(d[1] / l)
        L.append(l)
    Hxi = [[[-8, -4], [-4, 0]], [[0, 4], [4, 0]], [[0, -4], [-4, -8]]]
    Gd = _zeros(2, 3)
    for k in range(3):
        H = _matmul(Ji, _matmul(Hxi[k], _transpose(Ji)))
        bxx, bxy, byy = H[0][0], H[0][1], H[1][1]
        c, s = C[k], S[k]
        Gd[0][k] = Db_ * (c * bxx + (1 - nu) / 2 * c * byy + (1 + nu) / 2 * s * bxy)
        Gd[1][k] = Db_ * ((1 + nu) / 2 * c * bxy + (1 - nu) / 2 * s * bxx + s * byy)
    Aq = _zeros(3, 3)
    Bw = _zeros(3, 9)
    for k, (i, j) in enumerate(edges):
        for m in range(3):
            Aq[k][m] = L[k] * (C[k] * Gd[0][m] + S[k] * Gd[1][m]) / Ds
        Aq[k][k] -= 2.0 / 3.0 * L[k]
        Bw[k][3 * j] += 1.0
        Bw[k][3 * i] -= 1.0
        for n in (i, j):
            Bw[k][3 * n + 1] += L[k] / 2 * C[k]
            Bw[k][3 * n + 2] += L[k] / 2 * S[k]
    Q = _solve(Aq, Bw)
    Kb = _zeros(9, 9)
    for xi, eta in ((0.5, 0.0), (0.5, 0.5), (0.0, 0.5)):
        dPxi = [[4 * (1 - 2 * xi - eta), -4 * xi], [4 * eta, 4 * xi], [-4 * eta, 4 * (1 - xi - 2 * eta)]]
        dP = [[Ji[0][0] * d[0] + Ji[0][1] * d[1], Ji[1][0] * d[0] + Ji[1][1] * d[1]] for d in dPxi]
        Bb = _zeros(3, 9)
        for i in range(3):
            Bb[0][3 * i + 1] = dN[i][0]
            Bb[1][3 * i + 2] = dN[i][1]
            Bb[2][3 * i + 1] = dN[i][1]
            Bb[2][3 * i + 2] = dN[i][0]
        for k in range(3):
            row = [C[k] * dP[k][0], S[k] * dP[k][1], C[k] * dP[k][1] + S[k] * dP[k][0]]
            for a in range(3):
                for d in range(9):
                    Bb[a][d] += row[a] * Q[k][d]
        Kg = _matmul(_transpose(Bb), _matmul(Db, Bb))
        for a in range(9):
            for b in range(9):
                Kb[a][b] += A / 3 * Kg[a][b]
    Bg = [[sum(Gd[a][k] * Q[k][d] for k in range(3)) / Ds for d in range(9)] for a in range(2)]
    for a in range(9):
        for b in range(9):
            Kb[a][b] += A * Ds * (Bg[0][a] * Bg[0][b] + Bg[1][a] * Bg[1][b])
    ib = [6 * i + c for i in range(3) for c in (2, 4, 3)]
    sg = [1, 1, -1] * 3
    for a in range(9):
        for b in range(9):
            K[ib[a]][ib[b]] += sg[a] * sg[b] * Kb[a][b]

    # drilling penalty on rz_i - (v,x - u,y)/2
    kd = drill * G * h * A / 3
    for i in range(3):
        bd = [0.0] * 18
        bd[6 * i + 5] = 1.0
        for j in range(3):
            bd[6 * j + 1] -= 0.5 * dN[j][0]
            bd[6 * j] += 0.5 * dN[j][1]
        for a in range(18):
            for b in range(18):
                K[a][b] += kd * bd[a] * bd[b]

    # to global
    T = _zeros(18, 18)
    for blk in range(6):
        for a in range(3):
            for b in range(3):
                T[3 * blk + a][3 * blk + b] = R[a][b]
    return _matmul(_transpose(T), _matmul(K, T))


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _element_stiffness(X, E, nu, h, drill=None):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for i, Xi in enumerate(X):
        ops.node(i + 1, *Xi)
    ops.section('ElasticMembranePlateSection', 1, E, nu, h, 0.0)
    args = ['-drill', drill] if drill is not None else []
    ops.element('ShellDKMT', 1, 1, 2, 3, 1, *args)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.system('FullGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    # one (unloaded) step numbers the dofs and sizes the system; printA then
    # re-forms the tangent, so the factorization does not leak into A
    ops.analyze(1)
    A = ops.printA('-ret')
    neq = int(round(math.sqrt(len(A))))
    dofs = [d for n in (1, 2, 3) for d in ops.nodeDOFs(n)]
    return [[A[dofs[a] * neq + dofs[b]] for b in range(18)] for a in range(18)]


def _plate_center_deflection(eleType, a, h, n, E=1.0e4, nu=0.3, q=1.0):
    """simply supported (hard) square plate under uniform pressure, n x n squares cut in triangles"""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    tag = lambda i, j: i * (n + 1) + j + 1
    for i in range(n + 1):
        for j in range(n + 1):
            ops.node(tag(i, j), a * i / n, a * j / n, 0.0)
    ops.section('ElasticMembranePlateSection', 1, E, nu, h, 0.0)
    ele = 1
    area = {}
    for i in range(n):
        for j in range(n):
            n00, n10, n11, n01 = tag(i, j), tag(i + 1, j), tag(i + 1, j + 1), tag(i, j + 1)
            if (i + j) % 2 == 0:
                tris = [(n00, n10, n11), (n00, n11, n01)]
            else:
                tris = [(n00, n10, n01), (n10, n11, n01)]
            for t in tris:
                ops.element(eleType, ele, *t, 1)
                ele += 1
                for nd in t:
                    area[nd] = area.get(nd, 0.0) + (a / n) ** 2 / 2 / 3
    for i in range(n + 1):
        for j in range(n + 1):
            fix = [1, 1, 0, 0, 0, 1]
            if i in (0, n):
                fix[2] = 1
                fix[3] = 1    # beta_y = -rx = 0 along x = const
            if j in (0, n):
                fix[2] = 1
                fix[4] = 1    # beta_x = ry = 0 along y = const
            ops.fix(tag(i, j), *fix)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for nd, w in area.items():
        ops.load(nd, 0.0, 0.0, -q * w, 0.0, 0.0, 0.0)
    ops.system('BandGeneral')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    D = E * h ** 3 / (12 * (1 - nu ** 2))
    return -ops.nodeDisp(tag(n // 2, n // 2), 3) * D / (q * a ** 4)


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

X3D = [[0.3, -0.2, 0.1], [2.1, 0.4, 0.7], [0.8, 1.9, 1.3]]


@pytest.mark.parametrize('h', [0.01, 0.3, 1.0])
def test_stiffness_matches_reference(h):
    E, nu = 2.0e4, 0.2
    K = _element_stiffness(X3D, E, nu, h)
    Kr = ref_stiffness(X3D, E, nu, h)
    scale = max(abs(v) for row in Kr for v in row)
    err = max(abs(K[a][b] - Kr[a][b]) for a in range(18) for b in range(18))
    assert err / scale < 1.0e-10


def test_drill_option():
    E, nu, h = 2.0e4, 0.2, 0.3
    K = _element_stiffness(X3D, E, nu, h, drill=1.0e-2)
    Kr = ref_stiffness(X3D, E, nu, h, drill=1.0e-2)
    scale = max(abs(v) for row in Kr for v in row)
    err = max(abs(K[a][b] - Kr[a][b]) for a in range(18) for b in range(18))
    assert err / scale < 1.0e-10


def test_rigid_body_modes():
    E, nu, h = 2.0e4, 0.2, 0.3
    K = _element_stiffness(X3D, E, nu, h, drill=1.0)
    scale = max(abs(v) for row in K for v in row)
    modes = []
    for d in range(3):
        u = [0.0] * 18
        for i in range(3):
            u[6 * i + d] = 1.0
        modes.append(u)
    for axis in ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]):
        u = [0.0] * 18
        for i in range(3):
            t = _cross(axis, X3D[i])
            for d in range(3):
                u[6 * i + d] = t[d]
                u[6 * i + 3 + d] = axis[d]
        modes.append(u)
    for u in modes:
        f = [_dot(K[a], u) for a in range(18)]
        assert max(abs(v) for v in f) / scale < 1.0e-10


def test_thin_plate_matches_dkgt_and_kirchhoff():
    a, h, n = 1.0, 1.0e-3, 8
    w_mt = _plate_center_deflection('ShellDKMT', a, h, n)
    w_gt = _plate_center_deflection('ShellDKGT', a, h, n)
    assert abs(w_mt - w_gt) / w_gt < 1.0e-4
    assert abs(w_mt - 0.00406235) / 0.00406235 < 0.02


def test_thick_plate_matches_mindlin():
    # Navier solution of the first order shear deformation theory, k = 5/6,
    # nu = 0.3, a/h = 10: w D / (q a^4) = 0.004273 (Reddy)
    w = _plate_center_deflection('ShellDKMT', 1.0, 0.1, 16)
    assert abs(w - 0.004273) / 0.004273 < 0.01


def test_self_weight():
    rho = 2.5
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for i, Xi in enumerate(X3D):
        ops.node(i + 1, *Xi)
        ops.fix(i + 1, 1, 1, 1, 0, 0, 0)  # an all-fixed model has no equations to solve
    ops.section('ElasticMembranePlateSection', 1, 2.0e4, 0.2, 0.3, rho)
    ops.element('ShellDKMT', 1, 1, 2, 3, 1)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.eleLoad('-ele', 1, '-type', '-selfWeight', 0.0, 0.0, -1.0)
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear')
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.reactions()
    rz = sum(ops.nodeReaction(i + 1, 3) for i in range(3))
    area = 0.5 * math.sqrt(_dot(*(2 * [_cross(_sub(X3D[1], X3D[0]), _sub(X3D[2], X3D[0]))])))
    # the section's getRho() is the mass per unit area, rho * h
    w = rho * 0.3 * area
    assert abs(rz - w) / w < 1.0e-12
