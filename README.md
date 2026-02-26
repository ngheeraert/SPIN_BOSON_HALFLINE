# Spontaneous emission of Schrödinger‑cat **microwave wavepackets** in a high‑impedance waveguide  
*(Fortran / TDVP multi‑coherent‑state “multipolaron” solver)*

This repository contains the **Fortran code used to generate the core time‑domain results** for *spontaneous emission into a 1D superconducting waveguide* in the **ultrastrong‑coupling spin‑boson regime**, using a **time‑dependent variational principle (TDVP)** with a **multi‑coherent‑state (multipolaron) ansatz**. 

The central physical message (and what this code lets you compute) is:

> In a **high‑impedance superconducting transmission line**, where the effective light–matter coupling can reach order unity, spontaneous decay of a two‑level emitter does **not** radiate a narrowband “single photon”. Instead it produces a **short, broadband microwave burst** whose quantum state is well described as a **Schrödinger‑cat wavepacket** (a superposition of two macroscopically distinct coherent radiation envelopes). 

---

## 1) Microwave / circuit viewpoint (why the coupling becomes “large”)

In 3D free space, spontaneous emission is typically weak and perturbative. In superconducting circuits, you can engineer the **electromagnetic environment**: a 1D waveguide acts like a controlled continuum of microwave modes, and **Josephson‑junction chains** (or related metamaterials) can realize a **large characteristic impedance**.

A convenient dimensionless measure of coupling is  
\[
\alpha \sim \frac{Z}{R_K},
\]
where \(Z\) is the waveguide impedance and \(R_K = h/e^2\) the resistance quantum. Long JJ chains can reach \(Z \sim R_K\), hence \(\alpha \sim 1\). 

**Engineering translation:**  
- \(\alpha \ll 1\): the emitter sees a “gentle” 50 Ω‑like bath → narrowband emission, close to textbook lineshape.  
- \(\alpha \sim 1\): the line is so high‑impedance that the emitter’s decay becomes **non‑perturbative** → **strong line broadening**, multi‑photon content, and cat‑like radiation. 

---

## 2) Model solved by this code (what is being simulated)

### Hamiltonian (spin‑boson / waveguide‑TLS)
The code simulates a two‑level system (TLS) coupled to a bosonic continuum (waveguide modes). The Hamiltonian is written in the convention
- **TLS splitting** as a \(\sigma_x\) term, and
- **light–matter coupling** as a \(\sigma_z\) coupling to a field quadrature,

to emphasize the natural role of coherent states at strong coupling.

### Discretized waveguide (frequency bins)
The waveguide continuum is discretized into `nmode` modes (think: **frequency bins**). The code uses an Ohmic spectral density with exponential cutoff:
- mode frequencies: \( \omega_k = dk\,(k-\tfrac12),\; k=1..N\),
- spacing: \(dk = \omega_{\max}/N\),
- coupling envelope:
\[
g_k = \sqrt{2\,\alpha\,\omega_k\,dk\,e^{-\omega_k/\omega_c}}.
\]


A real‑space grid is also defined to compute field profiles via Fourier transforms (useful for time‑domain wavepackets, “impulse responses”, windowing, etc.). 

---

## 3) What “cat emission” means here (physics you can reproduce)

### (i) Broadband, time‑localized microwave radiation
At ultrastrong coupling, the emitted radiation becomes **spectrally broad** and correspondingly **localized in the time domain** (short wavepacket). 
This is the waveguide analogue of “a fast transient” in microwave engineering: large bandwidth ↔ short pulse.

### (ii) Multi‑photon content (not limited to one quantum)
The radiated field can contain **more than one photon on average** at \(\alpha \sim 1\), and the full many‑body state can involve an enormous effective Hilbert space—one motivation for the variational method. 

### (iii) Separation of time scales: fast energy relaxation vs slow decoherence
A key prediction is a large separation \(T_2 \gg T_1\): energy dumps into the line quickly, while qubit decoherence (and thus full “which‑path” information transfer to the wavepacket) can be parametrically slower.   
This leads to **partially coherent cats** at intermediate times.

### (iv) Phase‑space certification via Wigner tomography + filtering
Because the radiation is not monochromatic, one defines an **effective traveling mode** by applying an optimized temporal/spatial filter, removing the static near‑field “dressing cloud”, and then computing a **Wigner distribution** of that filtered mode.   
A cat shows two positive lobes (two classical coherent envelopes) and a negative interference region (“whiskers”). 

> **Signal‑processing analogy:** this is essentially *mode‑matched demodulation*: pick a window/function \(w(t)\) aligned with the outgoing pulse, project the field onto that mode, then analyze the resulting single‑mode state.

---

## 4) Numerical method (why it works in the ultrastrong regime)

### Multi‑coherent‑state (“multipolaron”) ansatz
The wavefunction is expanded as a *quantum superposition of multimode coherent states*:
```text
|Ψ(t)⟩ = Σ_{n=1..Ncs} [ p_n(t) |↑⟩ |{f_{n,k}(t)}⟩  +  q_n(t) |↓⟩ |{h_{n,k}(t)}⟩ ] .
```
This compactly represents large bosonic displacements and strongly non‑Gaussian radiation (cats), while keeping the variational dimension manageable. 

### TDVP equations + RK4 integration
TDVP gives a closed set of ODEs for the parameters \((p,q,f,h)\), which the code integrates with a 4th‑order Runge–Kutta scheme. 

### Convergence & scaling
- The TDVP error decreases rapidly with number of coherent states, typically \(\sim N_{cs}^{-2}\), and works across coupling regimes.   
- Internally, a key linear inversion step scales like \(N_{cs}^6\), so keeping \(N_{cs}\) modest is important. 

### Optional adaptive basis growth
The driver supports **adding coherent components on the fly** when an error indicator exceeds a threshold (useful when the wavefunction becomes more structured during emission). 

---

## 5) Repository contents (high level)

- `main.f90` — program entry point: parses CLI args, launches a trajectory.
- `system.f90` — model, mode discretization, TDVP equations, overlaps, observables.
- `output.f90` — RK4 evolution loop, adaptive basis logic, and data export.
- `lapackmodule.f90`, `inverse.f90` — LAPACK helpers for TDVP linear solves.
- `typedefs.f90`, `consts.f90` — kinds/constants.
- `bash.sh` — example cluster submission helper for a “quench” run. 

---

## 6) Build

### Dependencies
- A Fortran compiler (`gfortran`, `ifort`, `ifx`, …)
- BLAS + LAPACK (or vendor equivalents)

### Compile
```bash
make
```

This builds an executable named `mpol` (see `makefile`). If you compile manually, the README previously used:
```bash
gfortran -O3 -cpp -DDP \
  typedefs.f90 consts.f90 lapackmodule.f90 inverse.f90 system.f90 output.f90 main.f90 \
  -llapack -lblas -o mpol
```
(Use `-framework Accelerate` instead of `-llapack -lblas` on macOS if desired.)

---

## 7) Run

**Important:** the program exits if you provide **no** command‑line arguments.  
Create the output directory first:
```bash
mkdir -p data
```

### Common knobs
Key options (defaults shown in brackets) include:
- `-npini`  initial number of coherent states \(N_{cs}\) [1] 
- `-npadd`  number of additional coherent states to add adaptively [0] 
- `-al`     coupling \(\alpha\) [0.1] 
- `-del`    TLS splitting \(\Delta\) [0.1] 
- `-nm`     number of modes [400] 
- `-wc`     cutoff frequency \(\omega_c\) [1000] 
- `-wmax`   max discretized frequency \(\omega_{\max}\) [1] 
- `-tmax`, `-dt` final time and timestep  

### Example: adaptive basis growth (quick smoke test)
```bash
./mpol -npini 1 -npadd 6 -al 0.1 -del 0.1 -nm 400 -dt 0.05 -tmax 200 \
  -wc 1000 -wmax 1 -prep 2 -tref 0.2 -merr 1e-6 -p0 5e-6
```


### Example: cluster batch helper
`bash.sh` shows how to parameterize a “quench” and submit it to PBS/Torque. 

---

## 8) Outputs (what files you get)

All outputs are plain‑text `.d` files written to `data/` by default. 
The filename suffix `<param>` is an encoded parameter string produced by `parameterchar(sys)`.

### (A) Time‑series diagnostics
- `ErEN_<param>.d` : columns `(t, tdvp_error, energy, norm)` 
- `spinXYZ_<param>.d` : columns `(t, <σx>, <σy>, <σz>)` 
- `np_<param>.d` : columns `(t, Ncs)` (useful in adaptive runs) 

> Note: `ps_<param>.d` exists in the current snapshot but is reserved for per‑component weight diagnostics and is not populated consistently; for reliable weights, use the final‑state `ps_fst_<param>.d` below.

### (B) Final variational parameters (at `t = tmax`)
- `fks_fst_<param>.d` : for each mode \(k\), stores `ω(k)` followed by
  - real parts of `(f_n(k), h_n(k))` for all `n`, then
  - imaginary parts of `(f_n(k), h_n(k))` for all `n`. 

- `ps_fst_<param>.d` : stores the complex amplitudes `(p_n, q_n)` as
  - real parts for all `n`, then
  - imaginary parts for all `n`. 

### (C) Derived field observables
- `nx_fst_<param>.d` : photon density in real space (spatial profile of the traveling + bound contributions) 
- `nk_end_<param>_<wigxmin>_<wigxmax>.d` : a small table of \(n(k)\) values after optional spatial windowing (currently outputs first few modes). fileciteturn9file2L33-L35  

---

## 9) Wigner tomography / “mode‑matched” analysis (how to certify cats)

The paper’s phase‑space plots are obtained by:
1. **Filtering out** the static near‑field dressing cloud,
2. Defining an **effective traveling wavepacket mode** (temporal/spatial mode matching),
3. Computing the Wigner function \(W(x,p)\) of that effective mode. 

The codebase contains routines to do this (see `output.f90`: `wignerise`, `re_wignerise`, `printwigner`). In the current `main.f90`, these routines are not invoked by default; to reproduce Wigner plots you can add a call after the main trajectory (or create a small post‑processing driver that loads `fks_fst`/`ps_fst` and calls `re_wignerise`).

---

## 10) Quick plotting (Python)

Here is a minimal example to plot the emitted photon density \(n(x)\) and the qubit inversion \(\langle\sigma_z\rangle\):
```python
import numpy as np
import matplotlib.pyplot as plt

# replace <param> with the suffix printed by the program at runtime
nx = np.loadtxt("data/nx_fst_<param>.d")
spin = np.loadtxt("data/spinXYZ_<param>.d")

x, n_x = nx[:,0], nx[:,1]
t, sx, sy, sz = spin[:,0], spin[:,1], spin[:,2], spin[:,3]

plt.figure(); plt.plot(x, n_x); plt.xlabel("x"); plt.ylabel("n(x)"); plt.title("Photon density profile")
plt.figure(); plt.plot(t, sz); plt.xlabel("t"); plt.ylabel("<σz>"); plt.title("Qubit inversion")
plt.show()
```

---

## 11) How to cite

If you use or adapt this code for research results, please cite:

- N. Gheeraert, S. Bera, S. Florens (2017), *Spontaneous emission of Schrödinger cats in a waveguide at ultrastrong coupling*, **New J. Phys. 19, 023036**. 
