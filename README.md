# pdk2025

Neutrino interaction analysis framework for proton decay (PDK) sensitivity studies. Models detector response for the e+π⁰ and μ+π⁰ signal channels and neutrino backgrounds.

## Compilation

**Dependencies:** ROOT (≥ 6.x, with Physics, EG, Minuit2, ROOTDataFrame), EvtGen, Boost (program_options), C++23 compiler.

```bash
mkdir build && cd build
cmake ..
make            # build all targets
make si_check   # build a single target
```

The `EVTGEN_DATA_DIR` environment variable must point to the EvtGen data directory at cmake time, or be set beforehand so CMake can locate the EvtGen installation.

---

## Common Options

Most executables share the same command-line interface via `parse_command_line()`:

| Flag | Long form | Description |
|------|-----------|-------------|
| `-i` | `--input` | Input ROOT file(s) (repeatable / space-separated) |
| `-r` | `--input-corr` | Correction/friend ROOT file(s) |
| `-o` | `--output` | Output ROOT file (required) |
| `-g` | `--genie-mode` | Read GENIE `gRooTracker` tree instead of `outtree` |
| `-n` | `--no-fsi` | Disable final state interactions |
| `-m` | `--mu-pi-mode` | Analyse μ+π channel instead of e+π |
| `-t` | `--n-tagging` | Enable neutron tagging veto |
| `-e` | `--n-tagging-eff` | Neutron tagging efficiency (default: 0.25) |
| `-u` | `--unit` | Weight unit conversion factor (default: 1.0) |
| `-c` | `--external-capture-count` | Use external neutron capture count column |

---

## Executables

### `main`
Full neutrino background analysis pipeline for the e+π⁰ (or μ+π) channel. Reads a tracker event tree, runs detector simulation (smearing, ring merging, Michel counting), applies the signal selection, and produces histograms and a reduced TTree.

```bash
./main -i input.root -o output.root
./main -i file1.root file2.root -o output.root -m   # mu-pi mode
./main -i data.root -g -o output.root               # GENIE input
```

**Output:** histograms (ring multiplicity, invariant masses, angles, channel pie chart) + `outtree` TTree.

---

### `neutrinoevent`
Same pipeline as `main` but designed for weighted (reweighted) samples. Tracks weighted event counts at each cut stage and reports effective sample size (ESS).

```bash
./neutrinoevent -i input.root -o output.root
```

**Output:** histograms + TTree snapshot with weight columns.

---

### `neutrinoeventGenie`
Variant of `neutrinoevent` for native GENIE `gRooTracker` input trees.

```bash
./neutrinoeventGenie -i genie.root -o output.root
```

---

### `si_check`
Inspects single-pion FSI kinematics. Applies the full cut sequence and fills histograms of total momentum, total invariant mass, hadronic invariant mass W, and neutron capture count at each stage.

```bash
./si_check -i input.root -o output.root
./si_check -i input.root -o output.root -m   # mu-pi mode
```

**Output:** histograms `total_p`, `total_m`, `W`, `n_capture` (and `_limit` windowed variants) at each cut level.

---

### `eff_check`
Reconstruction efficiency analysis. Applies sequential topology cuts and records kinematic distributions (lepton energy, π⁰ energy, lepton–π⁰ angle, invariant masses) at each stage to evaluate per-cut efficiency.

```bash
./eff_check -i input.root -o output.root
```

---

### `capture_check`
Analyses neutron-capture-tagged events. Applies looser kinematic cuts (proton mass < 1.2 GeV, proton momentum < 1.2 GeV) and saves the passing events as a TTree for downstream analysis.

```bash
./capture_check -i input.root -o output.root
```

**Output:** `sample_event` TTree with columns `E_lepton`, `E_pi0`, `cos_theta_lepton_pi0`, `weight`, `n_capture`, `smear_proton_mass`, `smear_proton_momentum`, etc.

---

### `record_pdk`
Records reconstructed PDK signal events to a ROOT file. Optionally loads a correction/friend tree for initial-proton kinematics. Stores both truth and smeared versions of lepton and photon 4-momenta, plus kinematic fit results.

```bash
./record_pdk -i signal.root -o recorded.root
./record_pdk -i signal.root -r corr.root -o recorded.root
```

**Output:** `sample_event` TTree with `truth` and `smeared` `EventRec` branches, kinematic fit branches, and `total_weight`.

---

### `plot_from_record_common`
Post-processes output from `record_pdk`. Generates 3D histograms (E_lepton × E_π⁰ × cos θ) split by transparency (FSI-transparent vs. opaque) and proton momentum region.

```bash
./plot_from_record_common -i recorded.root -o plots.root
```

---

### `construct_th3d`
Builds 3D efficiency histograms of lepton energy × π⁰ energy × lepton–π⁰ opening angle from neutrino background samples.

```bash
./construct_th3d -i input.root -o th3d.root
```

**Output:** `TH3D` objects for full sample, FSI-transparent subset, and proton-momentum peak/tail regions.

---

### `slice_event`
Debug utility: prints kinematic details of the first 20 reconstructed events to stdout. No output file is written.

```bash
./slice_event -i input.root -o /dev/null
```

---

### `dump`
Debug utility: filters to 3-ring events with valid kinematic fits and prints full 4-momenta for the first 20 events.

```bash
./dump -i input.root -o /dev/null
```

---

### `genie_rate`
Integrates GENIE cross sections against a neutrino flux to compute interaction rates on hydrogen and oxygen targets.

```bash
./genie_rate -f flux.txt -g genie_xsec.root
./genie_rate -f flux.txt -g genie_xsec.root -i 0.2 -a 10.0
```

| Flag | Description | Default |
|------|-------------|---------|
| `-f` / `--flux-file` | Text file: energy (GeV) + 4 flux columns (νμ, ν̄μ, νe, ν̄e) | required |
| `-g` / `--genie-xsec-file` | GENIE ROOT cross section file | required |
| `-i` / `--emin` | Lower integration energy (GeV) | 0.1 |
| `-a` / `--emax` | Upper integration energy (GeV) | 25.0 |

**Output:** Console — integrated rates, CC/NC breakdown, average cross sections per target.

---

### `freeproton`
Generates isotropic free-proton decay events (p → l + π⁰) in the proton rest frame and writes them as a tracker-format TTree.

```bash
./freeproton [lepton_pdg] [n_events] [output.root]
./freeproton -11 1000000 signal_epi.root    # e+π⁰, 1M events
./freeproton -13 500000  signal_mupi.root   # μ+π⁰, 500k events
```

All arguments are positional and optional (defaults: PDG -11, 1 000 000 events, auto-named output).

**Output:** `outtree` TTree with branches `nparticles`, `P[5][4]`, `pdg[5]`, `status[5]`.

---

### `fluxint`
Integrates a flux file and prints the total flux to stdout.

```bash
./fluxint flux.txt
```

Input is a two-column text file (energy in GeV, flux value).

---

## Typical Workflow

```bash
# 1. Generate signal events
./freeproton -11 1000000 signal.root

# 2. Record reconstructed signal with kinematic fits
./record_pdk -i signal.root -o recorded_signal.root

# 3. Run neutrino background analysis
./neutrinoevent -i background.root -o background_out.root

# 4. Check FSI / single-pion systematics
./si_check -i background.root -o si_plots.root

# 5. Build 3D efficiency matrix
./construct_th3d -i background.root -o eff_th3d.root

# 6. Generate 3D signal plots for comparison
./plot_from_record_common -i recorded_signal.root -o signal_plots.root
```
