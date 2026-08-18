residual-collision-detection

""ci" (https://github.com/andrealo20/residual-collision-detection/actions/workflows/ci.yml/badge.svg)" (https://github.com/andrealo20/residual-collision-detection/actions/workflows/ci.yml)
""license: MIT" (https://img.shields.io/badge/license-MIT-blue.svg)" (LICENSE)
""language: C99" (https://img.shields.io/badge/language-C99-00599C.svg)" (https://en.wikipedia.org/wiki/C99)
""tests: 45" (https://img.shields.io/badge/tests-45%20passing-brightgreen.svg)" (tests/)
""sanitizers: ASan + UBSan" (https://img.shields.io/badge/sanitizers-ASan%20%2B%20UBSan-brightgreen.svg)" (.github/workflows/ci.yml)
""verified against Pinocchio" (https://img.shields.io/badge/verified%20against-Pinocchio-c2410c.svg)" (#verified-against-pinocchio)

Sensorless collision detection for robot arms, from the torque residual.

The library prefix is "qcollide_", kept short to avoid symbol collisions; the repository is named for what it does.

If the dynamic model of an arm is known, the torque required to produce a given motion can be predicted. The difference between that prediction and the measured torque contains information about external forces acting on the robot. A residual observer turns that difference into a collision signal without requiring a force sensor, artificial skin, or additional contact hardware.

"Residual following a contact" (docs/images/residual_trace.png)

The residual tracks the external torque. With a perfect model it settles on the applied value and returns to zero; with model error and realistic sensor effects it rides on a residual floor that any practical threshold has to clear first.

The result in one figure

"Detection sensitivity against realism" (docs/images/tradeoff.png)

In the representative model-error and sensor scenario evaluated here, the minimum detectable torque increases from 0.1 Nm to 4.8 Nm.

That loss of sensitivity is the main result of this project. It illustrates why practical collision detection depends not only on the observer itself, but also on model identification, friction compensation and the quality of the measured state.

Every newton metre of unmodelled torque contributes to the residual floor and reduces the margin available for detecting contact.

---

Why the naive approach does not work

The obvious formulation estimates external torque directly from the equation of motion:

τ_ext = τ_meas - M(q) q̈ - C(q, q̇) q̇ - g(q)

where:

- "τ_ext" is the external joint torque,
- "τ_meas" is the measured/applied joint torque,
- "M(q)" is the mass matrix,
- "C(q, q̇)" is the Coriolis matrix,
- "g(q)" is the gravity vector.

The problem is the acceleration "q̈".

Joint acceleration is normally not measured directly. Estimating it from position requires differentiating the encoder signal twice, and each differentiation amplifies high-frequency noise. At high sampling rates, even small position quantisation or encoder noise can therefore produce acceleration estimates whose error is larger than the contact torque being detected.

The momentum observer avoids "q̈" entirely.

Define the generalized momentum as:

p = M(q) q̇

Its derivative can be written as:

ṗ = τ + C(q, q̇)ᵀ q̇ - g(q) + τ_ext

There is no acceleration term.

This identity follows from the skew-symmetry property:

Ṁ - 2C  is skew-symmetric

which is why this repository is careful about which Coriolis matrix it computes.

Integrating the momentum equation and comparing the predicted momentum with the measured momentum gives a residual that converges toward the external joint torque using position, velocity and torque measurements only.

This is the momentum-observer formulation implemented in M1.

---

M0: the plant

"src/dynamics.c" implements an n-link planar arm: mass matrix, Coriolis matrix, gravity, forward and inverse dynamics, and RK4 integration.

A planar model keeps the algebra small enough to inspect independently, which matters because every later experiment depends on the correctness of this plant.

Verified against Pinocchio

<a name="verified-against-pinocchio"></a>

The plant is the reference against which the observer is evaluated, so agreement with self-written tests alone is not sufficient. A shared modelling error could make both the implementation and its tests agree while still representing the wrong dynamics.

"Pinocchio" (https://github.com/stack-of-tasks/pinocchio) provides an independent rigid-body dynamics implementation. The model is therefore cross-checked against it, in the same spirit that "qdsp" (https://github.com/andrealo20/qdsp) is checked against CMSIS-DSP.

pip install numpy pin
python3 tools/check_dynamics.py

Maximum absolute disagreement over 60 random configurations at 1, 2, 3 and 4 links:

Quantity| Max difference| Limit
Mass matrix "M(q)"| 5.3e-15| 1e-10
Gravity "g(q)"| 2.8e-14| 1e-10
Coriolis term "C(q, q̇) q̇"| 8.1e-09| 1e-05
Forward dynamics "q̈"| 5.6e-08| 1e-04

The tolerances differ for specific numerical reasons.

"M" and "g" are closed-form in both implementations and agree to machine precision. The Coriolis calculation in this repository uses central differences, which limits its numerical accuracy. Forward dynamics then propagates that error through the inverse mass matrix:

q̈ = M(q)⁻¹ [τ - C(q, q̇) q̇ - g(q)]

The measured amplification is 2.4× median at "n = 2" and 3.4× at "n = 4", with the condition number of "M" reaching 760. Its tolerance therefore cannot reasonably be tighter than the Coriolis calculation from which it is derived.

An earlier version used a "1e-8" forward-dynamics limit and failed for numerical rather than physical reasons.

The property the whole project rests on

There are infinitely many matrices "C" for which:

C(q, q̇) q̇

produces the correct Coriolis/centrifugal vector.

The momentum observer additionally requires a representation satisfying the rigid-body identity:

Ṁ - 2C  is skew-symmetric

Choosing a different valid representation of the Coriolis vector can therefore leave the ordinary dynamics correct while breaking the momentum observer.

The property is tested directly by numerically differentiating "M" along the trajectory.

Measurement| Result
Skew-symmetry residual of "Ṁ - 2C"| 4.1e-09 (matrix scale 6.1)
Mass matrix asymmetry| 0
Smallest Cholesky pivot over 200 poses| 0.032
Inverse/forward dynamics round trip| 6.4e-14

Integration quality

Energy conservation is checked by convergence order rather than against a single chosen threshold.

That distinction matters: a model containing a systematic error can still show small drift when the timestep is sufficiently small.

Step| Relative energy drift over 20 s
2.0e-3| 2.33e-05
1.0e-3| 2.51e-06
5.0e-4| 3.89e-07

Each halving reduces the measured drift by roughly an order of magnitude, consistent with high-order convergence. A separate trajectory-convergence test gives an error ratio of 18–20 against the ideal RK4 ratio of 16.

---

M1: the observer

"src/observer.c" implements the momentum residual:

r = K [ p - p(0) - ∫₀ᵗ (τ + C(q,q̇)ᵀ q̇ - g(q) + r) dt ]

For the ideal continuous-time observer, each joint satisfies:

ṙ = K (τ_ext - r)

For a constant external torque "T", the ideal response is:

r(t) = T [1 - exp(-K t)]

The continuous-time observer therefore has the time constant:

τ_observer = 1 / K

This does not mean that collision-detection latency is bounded below by "1/K". A threshold below the final residual can be crossed well before one time constant has elapsed.

The "1/K" interpretation is also a continuous-time result. In the discrete implementation, the sample period must remain sufficiently small relative to the observer dynamics; in other words, "K × dt" must be small enough for the chosen discretisation to represent the intended response accurately.

Everything required by the observer is measured or estimated from measurable quantities — position, velocity and torque — with no joint acceleration term.

Measured on a two-link arm at 10 kHz, with the observer given the same model as the plant:

Measurement| Result
Residual after 2 s with no contact| 0.0
Residual against applied 4.0 / −2.0 Nm| 4.0019 / −2.0004
Time constant at "K = 20 / 50 / 100"| 0.0498 / 0.0198 / 0.0098 s
Expected continuous-time "1/K"| 0.0500 / 0.0200 / 0.0100 s
Contact on joint 0: loaded / unloaded joint| 5.9955 / 0.0005 Nm
Worst residual during free swing| 0.0000 Nm

Two things that would otherwise fail silently

The Coriolis transpose.

The momentum equation requires:

C(q, q̇)ᵀ q̇

not:

C(q, q̇) q̇

Using the wrong one produces a detector that can appear correct while the arm is stationary but develops incorrect residuals during motion.

The test suite therefore exercises free motion specifically to catch this failure. Deliberately substituting "C" for "Cᵀ" increases the free-swing residual from 0.0000 to 7.2 Nm and causes three tests to fail.

That deliberate fault injection verifies that the tests actually detect the error they were designed to guard against.

The integration rule.

The integrand tracks "ṗ", which grows with joint velocity. Forward Euler therefore introduces an error proportional to both timestep and motion.

On a two-link arm driven to 19 rad/s, forward Euler produced a spurious 0.33 Nm residual on a joint with no external torque — already within the range of a light contact.

Using the previous and current available integrand values in a trapezoidal update substantially reduces the discretisation error:

Timestep| Spurious residual, forward Euler| Trapezoidal update
1.0e-3| 0.333 Nm| 0.022 Nm
5.0e-4| 0.165 Nm| 0.010 Nm
2.5e-4| 0.082 Nm| 0.005 Nm

The complete observer remains first order overall because the residual feedback is explicit. The trapezoidal integration therefore reduces the error constant rather than changing the global convergence order.

---

M2: from residual to decision

"src/detector.c" turns the residual into a collision decision.

A bare threshold is usually not enough: a single sample above it may be noise or a transient rather than evidence of sustained contact. The detector therefore requires the residual to remain above the trigger threshold for a configurable number of consecutive samples.

Once triggered, it releases only when the residual falls below a lower threshold.

Both mechanisms improve robustness at the cost of latency.

For an ideal first-order response:

r(t) = T [1 - exp(-K t)]

the time required to cross a threshold "r_th", provided "0 < r_th < T", is:

t_cross = -(1/K) ln(1 - r_th/T)

Detection latency therefore depends on:

- observer gain "K",
- contact torque "T",
- detection threshold "r_th",
- sampling period,
- number of required confirmation samples.

It is not bounded below by the observer time constant "1/K".

Measured on the two-link arm, perfect model, "K = 50" ("1/K = 20 ms"):

Threshold| Detection latency
0.5 Nm| 2.2 ms
1.0 Nm| 4.6 ms
2.0 Nm| 10.3 ms
4.0 Nm| 32.2 ms

The 1 Nm case illustrates the distinction directly: detection occurs in 4.6 ms even though the observer time constant is 20 ms.

Observer gain| Continuous-time "1/K"| Latency at 1 Nm
"K = 20"| 50 ms| 11.3 ms
"K = 50"| 20 ms| 4.6 ms
"K = 100"| 10 ms| 2.4 ms

Two further results are worth stating.

Minimum detectable torque: 0.9989 Nm at a 1.0 Nm threshold.

This value was found by bisection under the ideal-model experiment. For a sustained constant contact, the ideal residual converges to the applied external torque, so the asymptotic minimum detectable contact approaches the threshold.

Zero false alarms at every threshold tested, down to 0.05 Nm.

This is a baseline rather than a general claim about the method. With a perfect model and ideal measurements there is essentially no residual floor to trigger on. M3 deliberately removes that assumption by introducing friction and parameter error.

Short contacts: a measurable limitation

A sufficiently short contact can end before the residual reaches a chosen threshold.

For a rectangular contact of magnitude "T" and duration "d", the largest ideal residual reached before the contact ends is:

r_peak = T [1 - exp(-K d)]

Detection is therefore impossible for that event if:

r_peak < r_th

before accounting for any additional confirmation-sample requirement.

Measured at "K = 50", for which "1/K = 20 ms":

Contact duration| Duration / (1/K)| Outcome in this test
2 ms| 0.10| missed
5 ms| 0.25| detected
20 ms| 1.00| detected
100 ms| 5.00| detected

The 2 ms result does not mean that every 2 ms contact is fundamentally undetectable or that it would be missed at every possible threshold.

It means that, for the contact magnitude, observer gain, threshold, sampling rate and confirmation rule used in this experiment, the residual did not rise far enough for long enough to trigger the detector.

Increasing "K" makes the residual respond faster, but also makes the observer more sensitive to model and measurement errors. M3 and M4 quantify the practical consequence of that trade-off.

---

M3: what model error costs

Everything above used the same dynamic parameters in the plant and observer.

That is an idealised condition.

M3 introduces two controlled sources of modelling error that are hidden from the observer:

- friction, present in the plant only;
- inertial parameter error, implemented by scaling the observer's masses away from the plant values.

The residual floor

Even with no external contact, model mismatch produces a non-zero residual.

Any practical collision threshold must therefore sit high enough to avoid nuisance triggers caused by this residual floor.

Model error| Peak residual, no contact
none| 0.001 Nm
friction × 0.25| 0.603 Nm
friction × 0.50| 1.101 Nm
friction × 1.00| 1.866 Nm
inertia +2%| 0.178 Nm
inertia +5%| 0.445 Nm
inertia +10%| 0.890 Nm

In this experiment, parameter error enters the residual almost linearly: increasing the mass error from 5% to 10% raises the measured floor from 0.445 to 0.890 Nm, a ratio of 2.00.

The trade-off

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/bench/qcollide_sweep

For each threshold the benchmark asks two questions:

1. Does the detector fire when no external contact exists?
2. If not, what is the smallest tested contact it can reliably detect?

Case| Lowest usable threshold| Smallest contact caught there
perfect model| 0.10 Nm| 0.10 Nm
inertia 5% heavy| 1.00 Nm| 1.07 Nm
friction, representative case| 2.00 Nm| 3.61 Nm
both| 4.00 Nm| 5.28 Nm

Under the representative model errors evaluated here, the minimum detectable torque increases by roughly a factor of fifty relative to the ideal case.

A detector that catches a 0.1 Nm contact with a perfect model requires approximately 5.3 Nm in the combined friction and 5% mass-error experiment.

This experiment still excludes sensor noise, backlash and joint elasticity.

The result should therefore be interpreted as a controlled sensitivity study, not as a prediction of the detection threshold of a particular physical robot.

Two measurement mistakes worth recording

Reporting 0.00 Nm as excellent sensitivity.

When the threshold is below the residual floor, the detector is already firing without external contact. A bisection search for the minimum detectable torque can then converge toward zero.

An earlier implementation printed "0.00 Nm", which looked like excellent sensitivity while actually indicating an unusable threshold.

Those cases now report:

unusable (fires without contact)

Measuring at one instant.

The minimum detectable torque was initially measured with contact always starting at the same point in the trajectory.

Under friction, the answer depends strongly on contact phase. A collision arriving while the friction-induced residual is already close to the threshold requires less additional external torque to trigger the detector.

The benchmark now evaluates five contact phases and reports the worst case. This roughly doubled the measured value and makes the result substantially more representative of the tested trajectory.

---

M4: what the sensors cost

M3 gives the observer an imperfect model.

M4 instead introduces imperfect measurements, addressing the other major gap between an ideal simulation and implementation on a physical robot.

The momentum observer avoids acceleration estimation, but it still requires velocity.

On many robot joints, velocity itself is derived by differentiating encoder position.

For a 17-bit encoder sampled at 10 kHz, one encoder count corresponds to:

Δq = 2π / 2¹⁷ rad

Δq / Δt
= (2π / 2¹⁷) × 10⁴
≈ 0.48 rad/s

A one-count change between adjacent samples can therefore correspond to approximately 0.48 rad/s in the raw finite-difference velocity estimate.

That is large compared with the velocities encountered during slow motion.

"Sources of the residual floor" (docs/images/residual_floor.png)

Where the residual floor actually comes from

Measurement condition| Peak residual, no contact
exact velocity, no torque noise| 0.003 Nm
exact velocity, 0.05 Nm torque noise| 0.010 Nm
raw encoder difference| 26.313 Nm
encoder difference, 50 Hz filter| 0.477 Nm
filtered velocity + torque noise| 0.480 Nm

Two results stand out.

Raw differentiation of quantised encoder position makes the tested configuration unusable.

The resulting 26.313 Nm residual floor is much larger than the contacts of interest in these experiments.

Filtering reduces that floor by approximately a factor of 55.

In this experiment, torque noise contributes much less than velocity-estimation error.

Adding 0.05 Nm of torque noise on top of the filtered velocity estimate changes the measured residual floor from 0.477 to 0.480 Nm.

This result is specific to the sensor parameters and trajectory used here, but it shows why velocity estimation cannot be treated as a negligible implementation detail.

Filtering introduces its own trade-off.

A low-pass velocity estimate reduces quantisation noise but introduces lag. During acceleration, that lag appears to the observer as a modelling discrepancy.

Under an acceleration of approximately "4 rad/s²", the measured velocity estimate trails by 0.021 rad/s, compared with approximately 0.013 rad/s predicted from the filter time constant.

With everything present

Case| Lowest usable threshold| Smallest contact caught
perfect model, ideal sensors| 0.10 Nm| 0.10 Nm
real sensors, perfect model| 1.00 Nm| 0.88 Nm
friction + 5% parameters, ideal sensors| 4.00 Nm| 5.28 Nm
everything| 4.00 Nm| 4.80 Nm

These values belong to the specific simulated model, trajectory, friction parameters, sensor model and detection settings used by the benchmark.

They quantify how those effects interact; they are not claimed as universal limits for momentum-based collision detection.

A note on the roadmap

M4 was originally planned as a fixed-point implementation.

That direction was changed because the portfolio already investigates the cost of integer arithmetic in "qdsp" (https://github.com/andrealo20/qdsp) and "qcsi" (https://github.com/andrealo20/qcsi), while sensor realism addresses a question directly connected to this project's main objective.

A fixed-point observer remains a useful extension and is listed as open rather than silently removed.

---

What this does not do

These limitations are explicit because they materially affect how the results should be interpreted.

- Everything is simulated. The plant is a planar arm of uniform rods. No claim is made that the numerical thresholds measured here predict the performance of a particular physical robot.

- The friction parameters are representative, not identified from hardware. The Stribeck form is standard, but the coefficients used by the experiments were not fitted to a specific robot.

- The residual estimates external joint torque, not contact location. An n-joint robot produces n joint-torque residuals, whereas an arbitrary physical contact introduces both an unknown wrench and an unknown point of application. Collision detection is possible without uniquely solving that localisation problem.

- No backlash, joint elasticity or temperature drift is modelled. All can contribute additional residual error on real hardware.

- The sensor model is simplified. Encoder quantisation, finite-difference velocity estimation, low-pass filtering and torque noise are included, but the simulation is not intended to reproduce the complete sensing chain of a specific actuator.

- The first velocity sample is initialised from the simulated state. This avoids an artificial startup transient in the experiment. A hardware implementation using position encoders alone would instead need an explicit velocity-estimator initialisation strategy before enabling collision detection.

- No fixed-point implementation is provided. It was originally planned for M4 and replaced by sensor-realism experiments. The reasoning is documented in "docs/design.md".

---

Reproducing the numbers

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure

./build/bench/qcollide_sweep
./build/bench/qcollide_trace 0.0 0.00 0 > docs/data/trace_ideal.csv
./build/bench/qcollide_trace 1.0 0.05 1 > docs/data/trace_realistic.csv
python3 tools/make_figures.py

pip install numpy pin
python3 tools/check_dynamics.py

Every figure is generated from benchmark output rather than manually drawn, reducing the risk that plots and reported numerical results drift apart.

---

Build

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure

Requirements:

- CMake >= 3.16
- C99 compiler

Unity is fetched automatically.

---

Design notes

Implementation decisions, trade-offs and supporting measurements are documented in ""docs/design.md"" (docs/design.md), including a performance mistake that caused the sanitizer build to time out and how it was identified.

License

MIT — see "LICENSE" (LICENSE).