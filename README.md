# Cosmic Muon & Detector Response Simulator (C++ / CERN ROOT GUI)

An interactive graphical application built in **C++** and **CERN ROOT (`TGMainFrame`)** for simulating charged-particle energy loss, cosmic-ray muon fluxes at sea level, PMT threshold efficiencies, and 3D geometric acceptance in scintillation detectors and multi-layer cosmic telescopes.

## Key Features

### 1. Relativistic Bethe-Bloch & CSDA Range Calculator
* **Single-Medium & Inverse Multi-Layer Modes:** Computes linear stopping power ($-dE/dx$), mass stopping power ($-dE/\rho dx$), Bragg curves, and CSDA ranges integrated over 12,000 logarithmic energy steps. The inverse multi-layer solver (up to 10 layers in series $S_1 \dots S_N$) determines the minimum incident kinetic energy and momentum required for a particle to traverse the entire stack and stop at the bottom of $S_N$.
* **Physics Corrections:** Includes maximum kinematic energy transfer $T_{\max}$, Sternheimer–Peierls density effect parameterization, and Barkas–Berger shell corrections ($\eta = \beta\gamma \ge 0.1$).
* **Interactive Compound & Mixture Builder:** Dynamically constructs arbitrary molecules or mixtures (e.g., LAB scintillator, $\text{H}_2\text{O}$, PMMA) from base elements in `materials.txt`, evaluating molar mass, electron density, $\langle Z/A \rangle$, and Bragg-rule mean excitation energy $I_{\text{Bragg}}$.

### 2. Cosmic Muon Spectrum (Rastin 1984a) & PMT Threshold Efficiency
* Samples primary sea-level cosmic muons from the differential momentum spectrum of **Rastin (1984a)** with an angular distribution $I(\theta) = I_0 \cos^n(\theta)$.
* Propagates muons through passive overburden layers $S_1 \dots S_{N-1}$ into the active detector volume $S_N$ ($L_x \times L_y \times L_z$), separating **stopping muons** (candidates for $\mu \to e\nu\bar{\nu}$ decay measurements) from **through-going muons**.
* Applies a configurable photomultiplier (PMT) energy threshold ($\Delta E > E_{\text{thr}}$) to compute incident vs. detected rates and overall detection efficiency $\varepsilon_{\text{PMT}}$.

### 3. 3D Geometric Acceptance Monte Carlo & Plane Convergence Scan
* **`lab.cpp` (Single Parallelepiped Detector):** Generates $N_{\text{gen}}$ straight tracks across an ideal reference plane $S = L_x \times L_y$ at $z = Z$, classifying events into 5 topological categories (*Top→Bottom*, *Top→Side*, *Side→Bottom*, *Side→Side*, and *Miss*). Benchmarks MC rates, geometric factor $G_{\text{MC}}$, and mean chord length $\langle L \rangle$ against exact analytical formulas and generalized Cauchy–Crofton predictions.
* **`lab2.cpp` (3-Fold Coincidence Muon Telescope):** Models a 3-layer scintillator telescope geometry (Top active block + two intermediate/bottom trigger paddles separated by air gaps) to evaluate 3-fold coincidence rates and trigger efficiency.
* **Generation-Plane Convergence Scan:** Automatically sweeps 12 increasing generator plane sizes $S = L_{\text{gen}}^2$ to verify asymptotic plateau convergence $G_{\text{MC}}(S) \to G_{\text{Ana}}$.

## Repository Structure
* **`lab.cpp`** — Main GUI application featuring single-volume 3D geometric acceptance (5 topological classes), Bethe-Bloch/CSDA solver, compound builder, and Rastin/PMT spectrum simulator.
* **`lab2.cpp`** — Specialized GUI version configured for a 3-fold coincidence cosmic muon telescope geometry.
* **`materials.txt`** — Tab-separated database of elements, compounds, and scintillators (e.g., LAB, Polystyrene, Liquid Argon, Iron, Lead) with $Z/A$, mean excitation potential $I$, density $\rho$, and $Z_{\text{eff}}$.
* **`particles.txt`** — Tab-separated database of incident charged particles (muons, pions, kaons, protons, alphas, heavy ions).

## Prerequisites & Usage

### Requirements
* **CERN ROOT 6.x** compiled with GUI and 3D graphics support (`libGui`, `libGpad`, `libGraf3d`).
