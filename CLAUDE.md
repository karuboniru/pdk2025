# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

CMake with C++23. Standard out-of-source build:

```bash
mkdir build && cd build
cmake ..
make
# Or build a specific target:
make si_check
```

Dependencies: ROOT (Physics, EG, Minuit2, ROOTDataFrame), EvtGen, Boost (program_options). `EVTGEN_DATA_DIR` must be set in the environment.

The `DATA_PATH` macro is defined at compile time as `${CMAKE_SOURCE_DIR}/data`, pointing to the `data/` directory for pion angular/momentum distributions and EvtGen decay files.

There are no tests or lint commands.

## Architecture

### Shared Library (`event`)

All `src/*.cxx` files compile into the `libevent.so` shared library. The executables link against it.

**`NeutrinoEvent`** (`include/event.h`, `src/event.cxx`) — Core physics event class. Tracks particles at four stages:
- `in`: initial state (beam neutrino, target nucleon)
- `out`: primary interaction products
- `post`: post-FSI (final state interactions) particles
- `in_detector`: detector-level particles (as `momentum_pair` with `.truth` and `.smeared` fields)

Key method: `Rec_lpi_event(bool is_mu_pi)` reconstructs the lepton + pi0 (or pi+) from detector rings/showers. `finalize_and_decay_in_detector()` propagates post-FSI particles into detector rings and handles Michel electron counting.

**`FilterTrackedRDF`** (`include/commondefine.h`, `src/commondefine.cxx`) — Wraps ROOT's `RDataFrame::RNode` to track filter efficiency automatically. Use `FilterTracked()` instead of `Filter()` to record weighted event counts at each cut stage. Call `Report()` to print the cutflow table.

`DefineForEPi()` sets up derived columns for the e+π⁰ channel. `FilterSignalKinematics()` applies standard signal selection returning three signal region nodes.

**`EventRec`** (`include/data.h`) — ROOT-serializable struct (via `ClassDef`) for writing selected events to ROOT files. Contains lepton and gamma pair momenta plus kinematic fit results.

**Kinematic fitting** — Four variants in `src/kf*.cxx`: `kf` (base), `kf5dof`, `kf_egg`, `kf_full`. Results are stored in `EventRec`.

**`momentum_t`** — Alias for `ROOT::Math::PxPyPzEVector` (4-momentum). Used throughout.

### Executables

Each root-level `.cxx` file produces one executable:

| Executable | Purpose |
|---|---|
| `main` | Primary neutrino background analysis |
| `neutrinoevent` | Neutrino event reconstruction/analysis |
| `si_check` | Single-pion channel validation plots |
| `eff_check` | Detection efficiency analysis |
| `capture_check` | Neutron capture analysis |
| `record_pdk` | Record PDK signal events to ROOT file |
| `plot_from_record_common` | Post-processing plots from recorded events |
| `construct_th3d` | Build 3D efficiency histograms |
| `neutrinoeventGenie` | Analysis on GENIE-format input |
| `genie_rate` | GENIE interaction rate calculations |
| `slice_event` | Filter/slice events from ROOT files |
| `dump` | Debug dump of event contents |
| `freeproton` | Free proton interaction analysis |
| `fluxint` | Neutrino flux integration |

### Data Flow

Typical analysis pipeline:
1. Input ROOT/NTuple files contain neutrino interaction records
2. `FilterTrackedRDF` wraps an `RDataFrame` over the input
3. Per-event, a `NeutrinoEvent` is constructed and `finalize_and_decay_in_detector()` is called
4. Channel identification via `get_channelname()` / `get_channelname_no_nucleon()`
5. `DefineForEPi()` computes derived kinematic columns
6. `FilterSignalKinematics()` applies cuts; `Report()` prints cutflow
7. Histograms filled and saved to output ROOT files
