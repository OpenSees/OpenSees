# CFSTsteel material examples

CFSTsteel represents the cyclic uniaxial response of steel in concrete-filled steel tube members and bare steel members, including HSS and H-sections. It captures local buckling, combined hardening, and cyclic deterioration; calibrate the parameters for the section and confinement conditions. The material formulation is described by Wang, Wang and Lignos (2024), [DOI: 10.1002/eqe.4204](https://doi.org/10.1002/eqe.4204).

Both examples impose the same 1,200-step cyclic strain history and write strain, stress and tangent modulus to CSV. Use an OpenSees/OpenSeesPy build that includes this material:

```text
OpenSees CFSTsteel_cyclic.tcl
python CFSTsteel_cyclic.py
```

The examples use the original author parameter set, matching the symmetric animation in the documentation. These demonstrate the interface; they are not universal calibrations or reproductions of a specific experimental specimen. With `N=1`, only the first `Ck/gammak` pair is supplied. Strains are dimensionless and stress units are MPa. Output is written to the current working directory.

```text
uniaxialMaterial CFSTsteel tag E0 fy b2 sigma_cm a sigma_rm b3 \
    Lambda_t Lambda_cp Lambda_a Lambda_yp Lambda_r \
    c_t c_cp c_a c_yp c_r Qinf b N C1 gamma1 ... CN gammaN
```

Every scalar argument is required. `N` is an integer from 0 through 8, followed by exactly `N` pairs of kinematic hardening parameters. Compression is negative; use consistent units for `E0`, `fy`, `sigma_cm`, `sigma_rm`, `Qinf` and each `Ck`. The author's reference values `Lambda_t=11`, `Lambda_a=4`, and deterioration exponents equal to 1 are not inserted automatically by the parser.

The standalone C++ lifecycle regression uses the current OpenSees core classes and checks copying, rejected trials, reset, serialization and input errors. On a system with a C++17 compiler and LAPACK:

```text
cmake -S tests/cfststeel -B build-cfststeel
cmake --build build-cfststeel --config Release
ctest --test-dir build-cfststeel -C Release --output-on-failure
```

This model does not explicitly represent ductile crack initiation or propagation due to ultra-low-cycle fatigue. Softening response in fiber models also depends on the chosen element and integration discretization; reproduce the modeling assumptions of the cited study when comparing results.
