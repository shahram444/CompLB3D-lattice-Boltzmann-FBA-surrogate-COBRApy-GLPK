# CompLB3D

**CompLB3D: 3D pore-scale reactive transport with microbial metabolism**

A 3D pore-scale reactive transport model: lattice-Boltzmann flow and solute
transport coupled to microbial metabolism, aqueous speciation, and mineral
precipitation and dissolution that change the pore geometry as the simulation
runs.

CompLB3D extends [CompLaB](https://bitbucket.org/MeileLab/complab) from 2D to
3D and adds an equilibrium chemistry solver, an abiotic kinetics path, and
mineral dissolution.

---

## Authors

| | |
|---|---|
| **Shahram Asgari** | Department of Marine Sciences, University of Georgia, Athens, GA, USA — <shahram.asgari@uga.edu> |
| **Christof Meile** | Department of Marine Sciences, University of Georgia, Athens, GA, USA |

Developed in the Meile Lab, University of Georgia. CompLB3D builds on the
two-dimensional CompLaB, developed by Heewon Jung (Chungnam National
University) and co-workers with the University of Georgia.

---

## What it does

| | |
|---|---|
| **Flow** | Navier–Stokes by lattice Boltzmann (D3Q19), in an arbitrary pore geometry |
| **Transport** | advection–diffusion for any number of solutes (D3Q7) |
| **Biomass** | three interchangeable solvers: cellular automaton, finite difference, lattice Boltzmann |
| **Metabolism** | hand-written kinetics, flux balance analysis through **GLPK** or **COBRApy**, or a trained **surrogate** network — chosen per organism |
| **Speciation** | a components-and-species equilibrium tableau, solved in every pore voxel |
| **Geometry change** | mineral precipitation clogs pores; dissolution reopens them; the flow is re-solved either way |

Every one of those is optional and switched from `CompLaB.xml`. A run with all
of them off is a plain flow-and-transport solver.

---

## What changed in this version

The theme is that **the model no longer asks you to run anything before or
after it.** Preparing a metabolic model, building a pore space, fitting a
surrogate and checking mass balance were four separate scripts in Python and
MATLAB, run by hand, in the right order. They are now four blocks in
`CompLaB.xml`, handled natively in C++.

| | before | now |
|---|---|---|
| **Metabolic model** | download from BiGG, run `extractMM.py` to flatten it to a matrix file | `<model_source>bigg:iJO1366</model_source>`. Bundled, cached or downloaded, and checked against `models/manifest.txt`. `<model_filename>` also reads SBML directly now — the format is decided by the file's root element, so old matrix files still load |
| **Exchange reactions** | `<exchange_reaction_indices>27 35 -1</exchange_reaction_indices>` — column numbers, correct only for one revision of one model, and undetectable when wrong | `<exchange_reaction_names>EX_glc__D_e EX_o2_e none</exchange_reaction_names>` — resolved against the model at start-up. A wrong name stops the run and lists the near matches |
| **Pore space** | run `tools/geometry.py` to build a packing or threshold a volume | `<generate>spheres</generate>` inside `<domain>`, or `<import_raw>`. Inspected before the flow solver, written back so the run is reproducible from its own output, and a non-percolating domain stops the run |
| **Surrogate network** | `generateTrainingData.py` → `trainSurrogate.py` or `.m` → `exportSurrogateHeader.m` → paste the weights into `surrogateModel.hh` → **rebuild** | `<surrogate><train_if_missing>true</train_if_missing>`. The weights live in a file read at run time, so no rebuild; the training range travels with them and the solver clamps and counts instead of extrapolating |
| **Did it conserve?** | load the VTI volumes into Python afterwards and hope the answer is still in them | `<diagnostics>` writes `summary.csv` as the run goes and checks the sums you name, reporting drift the moment it appears |
| **Build** | edit the Palabos path inside `CMakeLists.txt` | `cmake -DPALABOS_ROOT=... ..`, or the environment variable. Stops with a clear message if it is not there |

None of the five blocks is required and every one is absent by default, so an
input file written before them takes exactly the path it always took. The
Python tools in `tools/` are still shipped and still supported: they are the
independent implementation the C++ was checked against, and they do things C++
cannot, such as reading a TIFF stack or drawing a plot.

`CompLaB_new_capabilities.xml` is a complete, runnable input file that turns all
of this on at once, with every parameter written out and commented.

### Verified against Palabos v2.3.0

This version has been **built and run**, which the previous one had not. Doing
that found four defects that no amount of reading would have:

| | |
|---|---|
| `EquilibriumChemistry<T>::MIN_CONC` | a `static constexpr` member with no out-of-class definition. It compiled cleanly and failed at the **link** step, because `std::min`/`std::max` take their arguments by reference |
| an out-of-bounds index in `complab.cpp` | **segfaulted any run with an LBM biofilm microbe**, minutes in, after the flow had converged. The counter advanced per LBM microbe over a vector holding only the planktonic ones |
| the generated example kinetics headers | omitted the `KineticsStats` namespace `complab.cpp` calls unconditionally, so **every example failed to build**. The structural validator passed them all, because it never compiled anything |
| the geometry generator | was handed the internal `nx`, which is two larger than the number in `<nx>`, and wrote a file two slices too wide. `readGeometry` stopped early and ignored the rest, leaving a truncated pore space that still looked plausible |

On `e_coli_core` the GLPK path reproduces the published COBRA values exactly —
0.8739 h⁻¹ aerobic on 10 mmol/gDW/h glucose, 0.2117 h⁻¹ anaerobic — and a
surrogate fitted to it during the run agrees with solving the linear program in
every voxel to four parts in ten thousand, about 150 times faster.

**Known, not fixed:** examples 13 to 15 report negative concentrations. With the
chemistry switched off entirely the same undershoot appears at iteration 0, so
it is the advection–diffusion Dirichlet inlet, in the transport core, and not
the reaction. `WHAT_TO_UPLOAD.md` has the measurements.

---

## Requirements

| | |
|---|---|
| Palabos | **v2.3.0**, downloaded separately (see below) |
| CMake | ≥ 3.10 |
| A C++11 compiler | GCC or Clang |
| MPI | optional but on by default |
| GLPK | optional, only for `-DENABLE_GLPK=ON` |
| Python 3 + cobrapy | optional, only for `-DENABLE_COBRAPY=ON` |

**Palabos is not in this repository.** It is a third-party library of about
100 MB under its own licence; vendoring it would make this repo unclonable on a
slow link and muddy the licence question. Download it from
<https://palabos.unige.ch/> and point CMake at it.

## Building

```bash
git clone https://github.com/shahram444/CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK.git
cd CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK

mkdir build && cd build

# Point CMake at your Palabos v2.3.0 tree. It can also come from the
# PALABOS_ROOT environment variable; without either, the default is the
# Tahoma path and the build stops with a message naming what it looked for.
cmake -DPALABOS_ROOT=/path/to/palabos-v2.3.0 ..                   # plain
cmake -DPALABOS_ROOT=/path/to/palabos-v2.3.0 -DENABLE_GLPK=ON ..  # with GLPK
cmake -DPALABOS_ROOT=/path/to/palabos-v2.3.0 -DENABLE_COBRAPY=ON ..
make -j
```

`CompLaB.xml` is read from the working directory, by that exact name.

```bash
cd ..
./complab              # or:  srun ./complab
```

`comp.sh` is a Slurm template; `examples/runExamples.slurm` submits the whole
example suite.

---

## Start here

```bash
python3 examples/makeExamples.py
./examples/runAllExamples.sh .
```

Sixteen cases, one per capability, each running in seconds on a laptop. Each is
a self-contained folder with its own `CompLaB.xml`, geometry, kinetics headers
and a README saying what to expect. `examples/README.md` is the map, and has a
table of which switches each metabolic case sets.

**Start with `examples/16_complete_pipeline`.** It is the shortest statement of
what this version does: a geometry and a `CompLaB.xml`, and from those it
unpacks a genome-scale model, resolves the exchange reactions by name, fits a
surrogate network to that model at start-up, and writes its own summary CSV.
Three seconds, no preparation step.

Four of them have answers that do not come from this code — an analytic
diffusion profile, mass balances, and a flux balance case whose optimum is
`min(Vs, 2·Vo)` by hand. Those are the ones to trust.

---

## Configuring a run

Everything is in `CompLaB.xml`. Two reference files document it:

| | |
|---|---|
| `CompLaB_reference_template.xml` | every tag the solver reads, with units, defaults and known limits |
| `CompLaB_new_capabilities.xml` | the four blocks added in this version, every parameter written out and explained. **Runnable as it stands** |

Two choices are per organism and independent of each other:

```xml
<microbe0>
    <solver_type>CA</solver_type>          <!-- how biomass moves: CA, FD, LBM -->
    <reaction_type>glpk</reaction_type>    <!-- how it metabolises -->
</microbe0>
```

`<reaction_type>` accepts `none`, `kinetics`, `glpk`, `cobrapy`, `surrogate`,
and the combinations `glpk_and_kinetics`, `cobrapy_and_kinetics`,
`surrogate_and_kinetics`, which run both paths and add the rates. Each also
needs its global switch on (`<enable_fba_glpk>` and friends), and a mismatch
stops the run with a message rather than doing something quietly wrong.

**Rate laws are compiled in.** `defineKinetics.hh` and
`defineAbioticKinetics.hh` are `#include`d, not read at run time, so changing
your chemistry means a rebuild. Every example carries its own pair.

---

## Fitting your own surrogate

**The short way is `<surrogate>` in the XML** — see the table above and
`CompLaB_new_capabilities.xml`. The run fits the network itself, from the
metabolic model, with no scripts and no rebuild.

The offline route is still here and is still the right one when you want to
compare architectures, train on a machine other than the one you will run on,
or keep the training data as a separate artefact. `surrogate_training/` has the
FBA sweep, a MATLAB trainer and a licence-free Python one that emit
byte-identical C++, and a verifier that compiles the generated header and checks
it against the fit.

`inspectSurrogate.py` reports any network's **valid input range**, recovered by
inverting the `mapminmax` scaling. Worth running before trusting a network
someone else fitted, including the one shipped in `surrogateModel.hh`: outside
its training box a network does not fail, it returns confident nonsense.

---

## Licence and third-party code

CompLB3D is **AGPL-3.0-or-later**, the same as CompLaB and Palabos. See
`LICENSE`.

`src/complab3d_glpkcpp.hh` derives from GLPKMEX / the Octave-Forge `glpkcc.cpp`
GLPK interface, which is GPL-licensed. It is included here on that basis and
its provenance is recorded in the file header. Anyone redistributing should
satisfy themselves that this is correct for their situation.

Palabos is a separate work under its own licence and is not distributed here.

---

## Citing

Use the **"Cite this repository"** button in the sidebar, which reads
`CITATION.cff`, or cite directly:

> Asgari, S. and Meile, C. (2026). *CompLB3D: a 3D pore-scale reactive
> transport model with microbial metabolism* (version 1.0.0).
> https://github.com/shahram444/CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK

```bibtex
@software{asgari_complb3d_2026,
  author  = {Asgari, Shahram and Meile, Christof},
  title   = {{CompLB3D: a 3D pore-scale reactive transport model with
             microbial metabolism}},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/shahram444/CompLB3D-lattice-Boltzmann-FBA-surrogate-COBRApy-GLPK}
}
```

Please cite the 2D CompLaB paper alongside this one: CompLB3D extends that work
rather than replacing it.

## Contributing

See `CONTRIBUTING.md`. The short version: run the example suite before opening
a pull request.
