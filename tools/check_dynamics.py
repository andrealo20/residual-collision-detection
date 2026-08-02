#!/usr/bin/env python3
"""
Cross-check the C dynamics against Pinocchio on identical inputs.

The plant in src/dynamics.c is the reference the collision observer is judged
against, so it cannot be taken on trust: an error there would show up as an
observer that appears to work while measuring the wrong thing. Pinocchio is
an independent implementation by another group, so agreement between the two
rules out a shared misunderstanding of the model, which is the failure a
self-written check cannot catch. This is the same arrangement qdsp has with
CMSIS-DSP.

The comparison builds the same planar arm twice: once as a Pinocchio model
assembled joint by joint, once by calling the C library through a small
harness. What is compared is the mass matrix, the gravity vector, the
Coriolis term and the forward dynamics, at random configurations.

Usage:
    python3 tools/check_dynamics.py

Requires: numpy, pin (Pinocchio). The C harness is built automatically.
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

HARNESS = r"""
#include "qcollide/dynamics.h"
#include <stdio.h>
#include <stdlib.h>

/* Reads: n mass length gravity q[0..n-1] qd[0..n-1]
   Writes: M row-major, then C row-major, then g, then qdd for zero torque. */
int main(int argc, char **argv) {
    qcollide_model m;
    double q[QCOLLIDE_MAX_DOF], qd[QCOLLIDE_MAX_DOF], tau[QCOLLIDE_MAX_DOF];
    double M[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double C[QCOLLIDE_MAX_DOF * QCOLLIDE_MAX_DOF];
    double g[QCOLLIDE_MAX_DOF], qdd[QCOLLIDE_MAX_DOF];
    size_t i, j, n;

    if (argc < 5) return 1;
    n = (size_t)atoi(argv[1]);
    if (qcollide_model_uniform(&m, n, atof(argv[2]), atof(argv[3]),
                               atof(argv[4])) != QCOLLIDE_OK) return 2;
    for (i = 0; i < n; ++i) q[i]  = atof(argv[5 + i]);
    for (i = 0; i < n; ++i) qd[i] = atof(argv[5 + n + i]);
    for (i = 0; i < n; ++i) tau[i] = 0.0;

    qcollide_mass_matrix(&m, q, M);
    qcollide_coriolis_matrix(&m, q, qd, C);
    qcollide_gravity(&m, q, g);
    qcollide_forward_dynamics(&m, q, qd, tau, NULL, qdd);

    for (i = 0; i < n; ++i) for (j = 0; j < n; ++j) printf("%.17g\n", M[i*n+j]);
    for (i = 0; i < n; ++i) for (j = 0; j < n; ++j) printf("%.17g\n", C[i*n+j]);
    for (i = 0; i < n; ++i) printf("%.17g\n", g[i]);
    for (i = 0; i < n; ++i) printf("%.17g\n", qdd[i]);
    return 0;
}
"""


def build_harness(tmpdir):
    src = os.path.join(tmpdir, "harness.c")
    exe = os.path.join(tmpdir, "harness")
    with open(src, "w") as f:
        f.write(HARNESS)
    cmd = ["gcc", "-std=c99", "-O2", f"-I{ROOT}/include", src,
           f"{ROOT}/src/dynamics.c", "-lm", "-o", exe]
    subprocess.run(cmd, check=True)
    return exe


def call_c(exe, n, mass, length, gravity, q, qd):
    args = [exe, str(n), repr(mass), repr(length), repr(gravity)]
    args += [repr(float(x)) for x in q] + [repr(float(x)) for x in qd]
    out = subprocess.run(args, capture_output=True, text=True, check=True)
    vals = [float(x) for x in out.stdout.split()]
    k = 0
    M = np.array(vals[k:k + n * n]).reshape(n, n); k += n * n
    C = np.array(vals[k:k + n * n]).reshape(n, n); k += n * n
    g = np.array(vals[k:k + n]); k += n
    qdd = np.array(vals[k:k + n])
    return M, C, g, qdd


def build_pinocchio(n, mass, length, gravity):
    """The same uniform-rod planar arm, assembled in Pinocchio.

    Links lie along local x, joints rotate about z, so the chain stays in the
    xy plane. Gravity is set along -y rather than the usual -z, because the
    C model is planar with gravity in the plane of motion.
    """
    import pinocchio as pin

    model = pin.Model()
    model.gravity.linear = np.array([0.0, -gravity, 0.0])

    parent = 0
    for i in range(n):
        placement = pin.SE3.Identity()
        if i > 0:
            placement.translation = np.array([length, 0.0, 0.0])
        joint = model.addJoint(parent, pin.JointModelRZ(), placement,
                               f"joint{i}")
        com = np.array([0.5 * length, 0.0, 0.0])
        # Uniform rod about its centre of mass, out-of-plane component
        # m*l^2/12; the in-plane components are irrelevant to planar motion
        # but must be positive for a valid spatial inertia.
        izz = mass * length * length / 12.0
        inertia_matrix = np.diag([1e-9, izz, izz])
        body = pin.Inertia(mass, com, inertia_matrix)
        model.appendBodyToJoint(joint, body, pin.SE3.Identity())
        parent = joint

    return model, model.createData()


def main():
    import pinocchio as pin

    rng = np.random.default_rng(0)
    mass, length, gravity = 1.3, 0.8, 9.81

    with tempfile.TemporaryDirectory() as tmp:
        exe = build_harness(tmp)

        overall = {"M": 0.0, "g": 0.0, "Cqd": 0.0, "qdd": 0.0}

        for n in (1, 2, 3, 4):
            model, data = build_pinocchio(n, mass, length, gravity)
            worst = {"M": 0.0, "g": 0.0, "Cqd": 0.0, "qdd": 0.0}

            for _ in range(60):
                q = rng.uniform(-np.pi, np.pi, n)
                qd = rng.uniform(-2.0, 2.0, n)

                M_c, C_c, g_c, qdd_c = call_c(exe, n, mass, length, gravity,
                                              q, qd)

                M_p = pin.crba(model, data, q)
                M_p = np.triu(M_p) + np.triu(M_p, 1).T   # crba fills upper only
                g_p = pin.computeGeneralizedGravity(model, data, q)
                # Pinocchio's nle is C(q,qd)qd + g(q); subtracting gravity
                # leaves the Coriolis term, which is what the C code returns
                # as a matrix and is compared here after multiplying by qd.
                nle_p = pin.nonLinearEffects(model, data, q, qd)
                Cqd_p = nle_p - g_p
                qdd_p = pin.aba(model, data, q, qd, np.zeros(n))

                worst["M"] = max(worst["M"], np.abs(M_c - M_p).max())
                worst["g"] = max(worst["g"], np.abs(g_c - g_p).max())
                worst["Cqd"] = max(worst["Cqd"], np.abs(C_c @ qd - Cqd_p).max())
                worst["qdd"] = max(worst["qdd"], np.abs(qdd_c - qdd_p).max())

            print(f"n = {n}:")
            for k in ("M", "g", "Cqd", "qdd"):
                print(f"  {k:4s} max |C - Pinocchio| = {worst[k]:.3e}")
                overall[k] = max(overall[k], worst[k])

        print()
        # Tolerances, and why they differ.
        #
        # M and g are closed-form in both implementations, so they agree to
        # machine precision and anything above 1e-10 is a real disagreement.
        #
        # The Coriolis term goes through central differences on M in the C
        # code, which caps its accuracy at around the square root of machine
        # epsilon times the scale of M. 1e-5 is loose enough for that and
        # tight enough to catch a wrong Christoffel form.
        #
        # qdd inherits the Coriolis error amplified by M^-1, so its tolerance
        # cannot be tighter than the Coriolis one. Measured amplification on
        # this arm: median 2.4 at n=2 and 3.4 at n=4, worst observed 24, with
        # condition numbers of M reaching 760. A limit of 1e-4 is the
        # Coriolis limit times an order of magnitude of headroom for that
        # amplification. An earlier version used 1e-8 here, which the code
        # failed for a reason that had nothing to do with correctness.
        limits = {"M": 1e-10, "g": 1e-10, "Cqd": 1e-5, "qdd": 1e-4}
        ok = True
        for k, lim in limits.items():
            status = "PASS" if overall[k] < lim else "FAIL"
            if overall[k] >= lim:
                ok = False
            print(f"  {k:4s} {overall[k]:.3e}  (limit {lim:.0e})  {status}")
        return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
