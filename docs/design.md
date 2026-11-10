# Design notes

## M0: The plant

### Why write the dynamics rather than call a library

Pinocchio would have supplied the model in a few lines. It is used here as a
*checker* instead, for a specific reason: the observer's correctness is
judged against this plant, so if the plant were a black box, a disagreement
between observer and plant would be untraceable. Writing it means every term
can be inspected, and using Pinocchio to verify it means a shared
misunderstanding of the model still gets caught.

Planar rather than spatial for the same reason. The algebra stays small
enough to check by hand, and nothing about the residual method depends on
the arm being planar.

### The Coriolis matrix, and why it is returned as a matrix

Most implementations return the Coriolis *vector* `C(q,qd) qd` directly,
which is cheaper and sufficient for simulation. The momentum observer needs
the matrix itself, because it uses the transpose: the identity

    pdot = tau + C(q,qd)' qd - g(q) + tau_ext

is what removes joint acceleration from the problem. Returning only the
vector would force the observer to reconstruct the matrix by finite
differences, defeating the purpose.

More importantly, *which* matrix is returned matters. Infinitely many C
satisfy `C qd = correct Coriolis vector`; only the Christoffel form makes
`Mdot - 2C` skew-symmetric, and the identity above depends on that. An
implementation that returned a different valid C would produce an observer
that runs, converges, and reports residuals that mean nothing. This is the
kind of error that has no symptom, so it is tested directly rather than
assumed: measured skew-symmetry residual 4.1e-09 against a matrix scale of
6.1.

### Finite differences for dM/dq, deliberately

The analytical derivative of the mass matrix is available but long, and an
error in it would break exactly the skew-symmetry above, silently. Central
differences are exact to second order, and the property they are meant to
produce is verified numerically. The cost buys a derivation that cannot be
quietly wrong.

### A performance mistake worth recording

The first version computed the derivatives inside the Christoffel triple
loop, calling a helper that evaluated the mass matrix twice per call. That
repeats the same handful of distinct derivatives O(n^3) times: for a 3-link
arm, 162 mass matrix evaluations per Coriolis call, and RK4 needs four of
those per integration step, 648 evaluations to advance the simulation by one
timestep.

There are only n distinct derivatives. Computing them once costs 2n
evaluations, twenty-seven times less work at n = 3, with the gap growing as
n^2.

How it was found is the part worth keeping: not by profiling, but because the
sanitizer build of the test suite ran long enough to hit a timeout. The
correctness tests were unaffected and would never have surfaced it. After the
fix the same suite runs in 22.8 seconds under ASan and UBSan, and a single
integration step costs 7.3 microseconds at n = 3.

### Tolerances derived rather than chosen

Every limit in `tools/check_dynamics.py` is justified in the file itself:

- M and g are closed-form on both sides, so 1e-10 flags a real disagreement.
- The Coriolis term is capped by finite-difference accuracy, hence 1e-5.
- `qdd` inherits the Coriolis error amplified by `M^-1`. Measured
  amplification: 2.4x median at n=2, 3.4x at n=4, worst observed 24, with
  cond(M) reaching 760. The limit is 1e-4.

That last one was originally 1e-8 and the check failed. The code was right;
the number was picked rather than derived. Measuring the amplification first
and setting the limit from it is the difference between a threshold that
means something and one that happens to pass.

### Energy: order, not magnitude

The obvious test asserts that energy drift stays below some value. That test
passes for a model with a sign error too, as long as the step is small
enough. Checking that drift shrinks with the step at the fourth-order rate
tests the integrator and the model together, and cannot be satisfied by
choosing a convenient step size.

## M1: The observer

### Why the residual feeds back into its own integral

The residual appears inside the integral that defines it. That looks
circular, and it is the point: it turns an open integration into a first-order
filter. An open integration of `tau + C' qd - g` would accumulate any model
error without bound, so a 1 % gravity error would produce a residual that
grows forever. With the feedback term the residual satisfies
`rdot = K (tau_ext - r)`, and a constant model error settles at a constant
offset instead of drifting. That is the difference between a detector that
needs periodic resetting and one that does not.

### The Coriolis transpose is the failure with no symptom

The momentum identity calls for `C'`, not `C`. The two are equal only when
the arm is at rest, so an implementation that uses the wrong one works on the
bench and fails in motion. Nothing crashes; the residual is simply wrong.

This is the error the test suite is built around. `test_residual_is_correct_while_moving`
drives the arm in free motion under gravity, where the Coriolis terms are
large and the true residual is exactly zero. To confirm the test actually
catches the mistake rather than merely passing, the transpose was removed on
purpose: the free-swing residual went from 0.0000 to 7.2 Nm and three tests
failed. A test that has never been seen to fail is not yet a test.

### Trapezoidal integration, found by a wrong number

The first version used forward Euler for the integral. The smoke test showed
the residual converging correctly on the loaded joint but sitting at 0.33 Nm
on a joint with no external torque, after the arm had been driven to 19 rad/s.

Halving the timestep halved the error, which identified it as first-order
integration error rather than a model mistake, the integrand tracks `pdot`,
and `pdot` grows with speed. The trapezoidal rule cuts it by fifteen times at
the same step.

Worth stating precisely: the scheme is still first order overall, because the
residual enters its own integrand explicitly and is therefore evaluated one
sample late. The trapezoidal rule improves the constant, not the rate. An
implicit treatment would recover second order at the cost of solving for r at
each step, which is not obviously worth it and has not been measured.

### What the observer does not do

The residual estimates external torque at the joints. Recovering the contact
force and where it acts is a different problem and is not attempted: an
n-joint arm yields n numbers, while a contact has six force components plus
an unknown location along the structure. Detection is well posed;
localisation needs assumptions this library does not make.

## M2: Detection

### Why the confirmation counter resets on a single sample

Requiring N consecutive samples above the threshold is the point of the
counter. An implementation that merely counted samples above the threshold
within a window would let uncorrelated noise accumulate towards a trigger
over an arbitrarily long time, turning a safeguard into a delay. The counter
therefore resets on the first sample below, and there is a test that runs
fifty near-misses in a row and asserts that nothing fires.

### Why release == threshold is rejected rather than allowed

It would be natural to treat equal thresholds as "hysteresis disabled". But a
residual sitting exactly at the threshold then produces a trigger and a
release on alternate samples, and a downstream safety system receiving a
1 kHz stream of contradictory decisions is worse off than with either
decision held consistently. The configuration is refused at init.

### Specifying contacts as joint torques, not forces at a point

A real contact is a force somewhere on a link, mapped to joint torques
through a Jacobian transpose. Skipping that mapping is deliberate: it makes
the ground truth exact, so a detection measurement reflects the detector
rather than the contact model. Adding a force-and-location interface would
add a second thing that could be wrong, and the residual method does not
depend on where the torque came from.

### The blind window is a property, not a bug

A rectangular contact of torque T and duration d drives the residual to at
most T(1 - exp(-Kd)), so a contact much shorter than the observer time
constant barely moves it. Measured at a 20 ms time constant and a 1 Nm
threshold: 2 ms missed, 5 ms detected.

This bounds the method, though not absolutely: a lower threshold or a larger
K recovers shorter contacts, and both pay by passing more model error and
noise, which is the trade-off M3 and M4 measure. Measuring it here means it
appears in the documentation as a stated limit rather than being discovered
by whoever deploys it.

### Zero false alarms is a baseline, not a result

With a perfect model, no friction and no sensor noise, the residual is
exactly zero when nothing touches the arm, so any threshold gives no false
alarms. Reporting that as a strength would be misleading. It exists to be the
control condition for M3, where friction makes the number non-zero and the
threshold finally has to be chosen against something.

## M3: Model error

### The asymmetry is the experiment

Friction is given to the plant and withheld from the observer. Inertial
parameters are perturbed in the observer's model and not the plant's. Neither
is a modelling convenience: that asymmetry is the entire experiment, because
it is what turns the threshold from a free parameter into a choice with a
cost on both sides.

Keeping the plant model and the observer's belief as separate objects, from
M1 onwards, is what made this possible without rewriting anything.

### Friction parameters are representative, not identified

The Stribeck form is standard in the identification literature. The values --
Coulomb around 0.8 Nm, stiction 1.2 Nm, viscous 0.5 Nm s/rad -- are of the
order reported for small collaborative-arm joints, but they were not
identified from a specific robot and the README says so.

That distinction matters for what can be claimed. The measurement shows how
sensitivity degrades as unmodelled torque of a realistic *shape and
magnitude* grows. It does not predict the numbers for any particular arm, and
presenting it as if it did would be the kind of overclaim this project exists
to avoid.

### Why the Coulomb term is regularised

The true Coulomb term is discontinuous at zero velocity. Left as a hard sign
function, an integrator stepping across zero sees an instantaneous torque
reversal, and RK4 responds by producing garbage rather than by failing. A
tanh over a narrow band agrees with the discontinuous version to machine
precision away from zero and removes the problem. The band width is a
parameter so its effect can be measured rather than assumed harmless.

### Driving the arm, not holding it

The sweep follows a sinusoidal trajectory. A detector evaluated at rest would
show almost no degradation, because friction and inertial error only appear
when the joints move, and it would report a trade-off that does not exist.

The first version commanded the arm with feedforward alone, inverse dynamics
of the sinusoid's acceleration, starting from rest. The sinusoid's velocity
is not zero at the start, so the arm began with a constant velocity error
that nothing corrected: without friction both joints wound up through
several radians in four seconds, and with friction the second joint stuck
and never moved. The friction rows were measured on a one-joint arm and the
frictionless rows on a spinning one, and the "phases" of the minimum
detectable torque were not phases of anything periodic. Every table looked
plausible. A review that simulated the motion on its own found it.

The arm now starts on the trajectory and tracks it with computed-torque
control from the plant model, at 40 rad/s bandwidth, acting on the measured
state. The gains do not enter the residual, since the observer is given the
torque actually commanded; they only decide how closely the arm follows the
path. Stiction still costs a few hundredths of a radian of tracking error at
each reversal. `tests/test_arm.c` pins the motion: tracking to 0.1 mrad without
friction, and both joints swinging at the commanded speed with it.

### Two mistakes in the measurement itself

**A zero that meant the opposite of what it read.** Below the residual floor
the detector already fires, so the bisection for minimum detectable torque
converges to zero, and the table printed "0.00 Nm" for exactly the thresholds
that were useless. It looked like the best row in the table. Those cases are
now reported as unusable, and the bisection is not run at all when the
false-alarm count is non-zero.

**A number that depended on luck.** The minimum detectable torque was first
measured with the contact always beginning at the same instant. Under
friction the residual is already oscillating, so a contact arriving in phase
with it needs far less torque to cross the threshold than one arriving out of
phase. The single-instant figure was the lucky case. Taking the worst over
five phases roughly doubled it -- from 0.15 Nm to 3.61 Nm in the
friction-only case -- and turned it into the number a detector actually has
to meet. Both of those figures were measured on the broken trajectory
described above. The sweep now tries ten phases and both directions of push,
holds each contact for 0.2 s rather than 1 s, and counts a detection only
while that contact can be blamed for it; the friction-only figure is 3.53 Nm.

Both errors flattered the result. That is the direction measurement mistakes
tend to run when the person measuring also wrote the thing being measured,
which is the argument for checking what a number would look like if the
system were broken.

## M4: Sensors

### Why this replaced the planned fixed-point milestone

The roadmap said M4 would implement the observer in fixed point. It was
changed deliberately. `qdsp` and `qcsi` already answer what integer
arithmetic costs, with a paired measurement in the second case, so a third
version would repeat a demonstration rather than add one. Sensor realism, on
the other hand, is directly on this project's thesis -- how well does
sensorless collision detection work, and what degrades it -- and it was the
last remaining difference between the simulation and a robot.

Changing a plan needs a reason rather than a preference, and the reason is
that M3 had already shown model error dominating, which made the question
"and what about the measurements" the obviously next one.

### The velocity estimate is the whole story

M1 avoided joint acceleration on the grounds that double differentiation
amplifies noise. The same argument applies with one differentiation less
force but it still applies, and it turns out to dominate everything else:

| Source | Residual floor |
|---|---|
| raw encoder difference | 27.840 Nm |
| filtered difference | 0.598 Nm |
| torque noise alone | 0.009 Nm |

Torque noise contributes under two percent of what velocity estimation does,
despite being the term that appears explicitly in the residual equation. The
intuitive ranking is wrong, and it would have stayed wrong if the two had
been introduced together instead of separately.

The filter is a trade rather than a fix. Measured lag under 4 rad/s^2 is
0.021 rad/s, against 0.013 predicted by the time constant alone -- the extra
comes from the difference operator itself, which also lags by half a sample.
On an accelerating joint that lag reaches the residual as a torque that is
not there, so a slower filter does not simply give a lower floor.

### The noise generator is explicit, not rand()

A deterministic xorshift with a caller-supplied seed rather than the C
library generator, because `rand()` differs between platforms and a sweep
whose numbers move between machines cannot be compared against itself. The
generator's moments are tested: a Box-Muller with the wrong constant produces
samples that look Gaussian and have the wrong variance, which would rescale
every noise figure in the project by a factor nobody would notice.

### A test that failed while the code was right

`test_noise_is_deterministic` compared two identically seeded generators with

    TEST_ASSERT_EQUAL_DOUBLE(gaussian(&a), gaussian(&b));

and failed. The generator was correct. Unity expands that macro with the
`expected` argument appearing twice -- once to compute the tolerance and once
for the comparison -- so `a` advanced twice per iteration and `b` once.

The rule this leaves: never call a stateful function inside an assertion
macro. Pull the values into locals first. The cost of getting it wrong is a
test that fails for a reason unrelated to what it is testing, which sends the
search into the wrong file entirely.

## Still open

- **A fixed-point observer.** Planned as M4 and replaced by sensor realism;
  see above. The scaling design is the interesting part, since the physical
  quantities here span a far wider dynamic range than the signals in qcsi.
- **Contacts on joint 1, and contacts on both joints at once.** The minimum
  detectable torque is measured with the contact on joint 0, in both
  directions. The floor differs per joint, so the other joint may be harder.
- Everything is simulated. No claim is made here about behaviour on a real
  arm, and the README repeats the limitation rather than leaving it implicit.
