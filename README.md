# Spontaneous emission of Schrodinger-cat wavepackets (Fortran / TDVP multi-coherent-state)

This repository contains the Fortran implementation used to generate the time-domain results for **spontaneous emission into a 1D bosonic waveguide / continuum** in the **ultrastrong-coupling spin–boson setting**, using a **time-dependent variational principle (TDVP)** with a **multi-coherent-state (multipolaron) ansatz**.

The code evolves a variational wavefunction of the form

```
|Psi(t)> = sum_{n=1..Np} [ p_n(t) |up>  |{f_{n,k}(t)}>  +  q_n(t) |down> |{h_{n,k}(t)}> ]
```

where `|{f_{n,k}}>` and `|{h_{n,k}}>` are multimode bosonic coherent states (one complex displacement per mode `k`).

The key idea is that increasing the number of coherent components `Np` allows the wavefunction to represent strongly non-Gaussian photonic states (including cat-like wavepackets), while TDVP yields a closed set of ODEs for the parameters `(p,q,f,h)`.

A reference describing the physics and the variational method is provided in the attached paper:

- *N. Gheeraert et al.* (2017), **“Spontaneous emission of Schrodinger cats in a waveguide at ultrastrong coupling”**.

---

## Contents

- `main.f90` – program entry point
- `system.f90` – model definition, mode discretisation, TDVP equations of motion, overlaps and cached mode sums, observables
- `output.f90` – time evolution driver (RK4), optional “add-polaron” growth, and output routines
- `typedefs.f90`, `consts.f90` – numerical kinds and constants
- `lapackmodule.f90` – thin wrappers around BLAS/LAPACK routines used for inversions and linear solves
- `inverse.f90` – helper for solving the TDVP “kappa” linear system by flattening a 3-index tensor

---

## Physical model (high-level)

The simulation corresponds to a two-level system (TLS) coupled to a bosonic continuum (waveguide modes). The discretised mode parameters are:

- Number of modes: `nmode`
- Frequencies: `w(k) = dk * (k - 1/2)`, `k = 1..nmode`
- Mode spacing: `dk = wmax / nmode`
- Ohmic coupling (discretised spectral density with exponential cutoff):

```
 g(k) = sqrt( 2 * alpha * w(k) * dk * exp( -w(k)/wc ) )
```

where:
- `alpha` controls the overall coupling strength,
- `wc` is the exponential cutoff frequency,
- `wmax` sets the maximum discretised frequency.

The code also defines a real-space grid used for Fourier transforms of the field:

- `length = pi / dk`
- `dx     = length / nmode`

These are used to compute field profiles and photon densities in position space.

---

## Numerical method (high-level)

### TDVP equations

The routine `CalcDerivatives(sys, st)` computes time-derivatives of all variational parameters:

- amplitudes: `pDot(:)`, `qDot(:)`
- coherent displacements: `fDot(:,:)`, `hDot(:,:)`

The implementation relies heavily on:

- **coherent-state overlap matrices** (`ov_ff`, `ov_hh`, `ov_fh`, `ov_hf`)
- cached mode sums (`bigW_*`, `bigL_*`) to avoid recomputing expensive sums over modes

A central step is solving a linear system for an auxiliary matrix (often called a “kappa” object in TDVP derivations). In this code, the relevant 3-index tensor is flattened into a `(Np^2 x Np^2)` system and solved with LAPACK.

### Time integration

Time propagation uses a classic **4th-order Runge–Kutta (RK4)** integrator (`evolve_RK4`).

The driver `evolveState_HL`:

- advances until `tmax`,
- periodically records observables into a `traj` structure,
- optionally **adds additional coherent states** ("polarons") when an error indicator exceeds a threshold and a minimum time has elapsed.

### Error monitoring / adaptive growth

The code computes a TDVP error indicator (`error(sys, ost, st)` in `system.f90`). In adaptive runs (`npadd > 0`), once:

- `t >= tref_arr(current_Np)` and
- `error >= merr_arr(current_add_index)`

the evolution stops, a new coherent component is appended (`addpolaron`), and evolution resumes.

---

## Building

### Requirements

- A Fortran compiler (e.g. `gfortran`, `ifort`, `ifx`)
- BLAS + LAPACK (or a vendor equivalent such as Apple Accelerate / Intel MKL)

### Compile (gfortran + system BLAS/LAPACK)

From the directory containing the `.f90` files:

```bash
gfortran -O3 -cpp -DDP \
  typedefs.f90 consts.f90 lapackmodule.f90 inverse.f90 system.f90 output.f90 main.f90 \
  -llapack -lblas -o cats_tdvp
```

Notes:
- `-cpp` enables preprocessing.
- `-DDP` selects double precision kinds (`typedefs.f90`). Remove it for single precision.

### Compile on macOS (Apple Accelerate)

```bash
gfortran -O3 -cpp -DDP \
  typedefs.f90 consts.f90 lapackmodule.f90 inverse.f90 system.f90 output.f90 main.f90 \
  -framework Accelerate -o cats_tdvp
```

---

## Running

**Important:** the program exits if you provide **no** command-line arguments. Pass at least one option.

### Common parameters

Key options (defaults shown in brackets):

- `-npini`  initial number of coherent states `Np` [1]
- `-npadd`  number of coherent states to add adaptively [0]
- `-al`     coupling strength `alpha` [0.1]
- `-del`    TLS splitting `Delta` [0.1]
- `-nm`     number of modes `nmode` [400]
- `-dt`     base time step [0.05]
- `-tmax`   final time [200]
- `-wc`     cutoff frequency [1000]
- `-wmax`   maximum discretised frequency [1]
- `-prep`   initial TLS state: `1 = (|up>+|down>)/sqrt(2)`, `2 = |up>`, `3 = |down>` [2]
- `-path`   output path selector: `d` -> `data/`, `.` -> current directory [`data/`]

Adaptive-growth knobs (only used when `-npadd > 0`):

- `-tref`   minimum time before adding a new coherent state [0.2]
- `-merr`   error threshold for adding a new coherent state [1e-6]
- `-p0`     initial amplitude assigned to the newly added component [5e-6]

### Example: fixed-basis run

```bash
./cats_tdvp -npini 1 -npadd 0 -al 0.1 -del 0.1 -nm 400 -dt 0.05 -tmax 200 -wc 1000 -wmax 1 -prep 2
```

### Example: adaptive basis growth

```bash
./cats_tdvp -npini 1 -npadd 6 -al 0.1 -del 0.1 -nm 400 -dt 0.05 -tmax 200 \
  -wc 1000 -wmax 1 -prep 2 -tref 0.2 -merr 1e-6 -p0 5e-6
```

---

## Outputs

All outputs are plain-text `.d` files written to `data/` by default.

### Time series

- `ErEN_<param>.d` : columns `(t, tdvp_error, energy, norm)`
- `spinXYZ_<param>.d` : columns `(t, <sigmaX>, <sigmaY>, <sigmaZ>)`
- `np_<param>.d` : columns `(t, Np)` (useful for adaptive runs)
- `ps_<param>.d` : columns `(t, |p_1|^2, |p_2|^2, ... )` (size matches `npini+npadd`)

### Final variational parameters (at `t = tmax`)

- `fks_fst_<param>.d` : for each mode `k`, stores `w(k)` followed by
  - real parts of `(f_n(k), h_n(k))` for all `n = 1..Np`, then
  - imaginary parts of `(f_n(k), h_n(k))` for all `n = 1..Np`.

- `ps_fst_<param>.d` : stores the complex amplitudes `(p_n, q_n)` as
  - real parts for all `n`, then
  - imaginary parts for all `n`.

### Derived field observables

- `nx_fst_<param>.d` : photon density in real space (computed from the up-component)
- `nk_end_<param>_<wigxmin>_<wigxmax>.d` : a small table of `n(k)` values (currently first 5 modes) after optional time-window filtering
- `fxs_*`, `fks_*` : Fourier-transformed or filtered field profiles (used internally for post-processing)

**Parameter strings:** `<param>` is a compact identifier produced by `parameterchar(sys)` that encodes the run parameters.

---

## Notes & tips

- **Performance scaling:** the TDVP linear algebra grows quickly with the number of coherent components `Np` (matrix/tensor sizes scale as `Np^2` and worse). Expect adaptive runs to become much heavier as `Np` increases.
- **Stability checks:** `checkTimestep` can attempt to “repair” steps when the instantaneous energy drift exceeds `max_deltaE`.
- **Precision:** compile with `-DDP` unless you have a strong reason not to.

---

## How to cite

If you use or adapt this code for research results, please cite the accompanying paper:

- N. Gheeraert et al. (2017), *Spontaneous emission of Schrodinger cats in a waveguide at ultrastrong coupling*.
