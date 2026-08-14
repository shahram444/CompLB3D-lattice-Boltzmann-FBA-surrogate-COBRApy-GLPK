# CompLB3D example suite

Sixteen tiny cases, one per capability. Each runs in seconds on a laptop, each
is self-contained, and each says what it should produce.

They are deliberately small. The point is to prove a feature works and to show
what the switches do, not to be physically interesting. Every domain is
24 × 24 × 6 with z periodic, which makes it quasi-2D: every z-plane sees
identical conditions, so a result that depends on z is a bug rather than
physics.

```bash
python3 makeExamples.py                       # writes all 16 case folders
python3 validateExamples.py ../CompLB3D       # structural check, one second
./runAllExamples.sh ../CompLB3D               # build and run everything
./runAllExamples.sh ../CompLB3D 09 16         # or just these two
```

The geometries are built by `makeExamples.py` itself. There used to be a
separate `makeGeometry.py` that wrote an `input_shared/` folder for it to copy
from; that was one more command to remember and nothing else ever read those
files.

---

## The cases

| | case | what it exercises |
|---|---|---|
| 01 | `flow_only` | Navier–Stokes plus advection–diffusion, no chemistry. The smoke test. |
| 02 | `diffusion_only` | `<Peclet>0`, so no flow. **Has an analytic answer.** |
| 03 | `abiotic_kinetics` | `defineAbioticKinetics.hh`: A + B → C |
| 04 | `equilibrium` | the speciation tableau, on the carbonate system |
| 05 | `biotic_cellular_automaton` | Monod growth, CA biomass solver |
| 06 | `biotic_finite_difference` | the same case, FD solver |
| 07 | `biotic_lattice_boltzmann` | the same case, LBM solver |
| 08 | `two_microbes` | two organisms, two *different* solvers, competing |
| 09 | `fba_glpk` | flux balance analysis through GLPK. **Has an analytic answer.** |
| 10 | `fba_cobrapy` | the same model through COBRApy. Must agree with 09. |
| 11 | `surrogate` | the trained network instead of an LP |
| 12 | `mixed_reaction_types` | `glpk_and_kinetics`: both paths, rates added |
| 13 | `precipitation` | mineral growth and pore clogging |
| 14 | `dissolution` | mineral loss and pore reopening |
| 15 | `precip_and_dissolution` | both directions on one mineral |
| 16 | `complete_pipeline` | `<model_source>`, `<exchange_reaction_names>`, in-run surrogate training, `<diagnostics>`. **Nothing prepared in advance.** |

Cases 04, 13, 14 and 15 also write `output/summary.csv` through `<diagnostics>`.
For 13, 14 and 15 the porosity column in that file *is* the result: pore space
sealing as the mineral grows and reopening as it dissolves.

Cases 05, 06 and 07 are **the same case with one line changed**. So are 09 and
10. Those pairs are the useful ones: they isolate a single decision.

---

## Choosing a metabolic solver

This is the question the sixteen files answer only by example, so here it is
once. **Two gates, and both have to agree.**

The global switch says *this solver is compiled in and available in this run*.
`<reaction_type>` inside each `<microbeN>` says *this organism uses it*. Asking
in one place and not the other stops the run with a message naming what is
missing; that is deliberate, so you cannot quietly get different biology from
the one you asked for.

| | 09 | 10 | 11 | 12 | 16 |
|---|---|---|---|---|---|
| `<enable_fba_glpk>` | true | false | false | true | false |
| `<enable_fba_cobrapy>` | false | true | false | false | false |
| `<enable_surrogate>` | false | false | true | false | true |
| `<enable_kinetics>` | false | false | false | true | false |
| `<reaction_type>` | `glpk` | `cobrapy` | `surrogate` | `glpk_and_kinetics` | `surrogate` |
| cmake | `-DENABLE_GLPK=ON` | `-DENABLE_COBRAPY=ON` | — | `-DENABLE_GLPK=ON` | `-DENABLE_GLPK=ON` |

`<reaction_type>` accepts `none`, `kinetics`, `glpk`, `surrogate`, `cobrapy`,
and the three combined forms `glpk_and_kinetics`, `surrogate_and_kinetics`,
`cobrapy_and_kinetics`, which run both paths and add the rates. The numeric
spellings 0 to 7 also work, in that order.

Case 16 needs GLPK even though no organism uses it, because it **trains** its
surrogate by sweeping a linear program at start-up. Once `output/ecoli.srg`
exists, the same case loads it and needs no solver at all.

### What each solver then demands of the input file

| | GLPK | COBRApy | surrogate |
|---|---|---|---|
| a metabolic model | yes | yes | only to train one |
| `<model_filename>` **or** `<model_source>` | yes | yes | only to train one |
| exchange mapping | `<exchange_reaction_names>` (SBML) or `<exchange_reaction_indices>` | same | same, when training |
| `<fba_maximum_uptake_flux>` | yes | yes | yes |
| `<substrate_lower_bounds>` / `_upper_bounds` | yes | yes | no |
| `<objective_direction>` | yes | yes | no |
| a fitted network | no | no | `<weights_file>`, or the weights compiled into `surrogateModel.hh` |

**Prefer `<exchange_reaction_names>` to `<exchange_reaction_indices>`.** An
index that points at the wrong reaction is undetectable: the linear program is
still well posed and the answer is wrong everywhere. A name that does not
resolve stops the run and prints the near matches. Names need an SBML model;
the `extractMM.py` matrix format does not carry them.

---

## What is worth checking, not just running

`runAllExamples.sh` checks that each case runs, exits cleanly and never
reports a negative concentration. That catches crashes and instabilities. It
does **not** check the numbers, and four of these cases have a right answer
that does not come from this code:

**02** — with no reaction and one face held at 1, the other at 0, the steady
profile is a **straight line**. Curvature means a reaction is firing somewhere
it should not, or a boundary is not being held.

**03, 13, 14, 15** — **mass balance**, with one distinction worth making
before you check it.

The RATE stoichiometry is exact by construction: one A plus one B makes exactly
one C, and one mole of calcite removed appears as exactly one mole of Ca²⁺.
That is what the rate law says and it is worth confirming.

The TOTAL of a species is a different claim, and it only holds for a species
nothing supplies from outside. In 03 both A and B enter at Dirichlet inlets, so
A + C rises for a reason that has nothing to do with the chemistry; the same is
true of Fe²⁺ in 13 and 15. **The one genuinely closed sum in the suite is
`Ca + calcite` in case 14** — calcite dissolves into Ca²⁺ and neither leaves the
domain.

These four now write `output/summary.csv` through `<diagnostics>`, so the totals
are a column rather than a post-processing job. To have the code check the
closed sum for you, add one line to case 14:

```xml
<diagnostics>
    ...
    <conserve>Ca+calcite</conserve>
</diagnostics>
```

That check is left commented rather than shipped on because it has not been run
here, and a check that fails on its first outing tells you nothing about your
input file. Turn it on and the run reports the drift the moment it exceeds
`<tolerance>`.

**09 and 10** — the toy model is four reactions, so the answer is a formula:

```
growth = min( Vs, 2 · Vo ),    V = Vmax · C / (Kc + C)
```

With Vmax 10 and 4 and both concentrations well above their half-saturation
constants, the acceptor limits and growth should settle near **8 mmol/gDW/h**.
That number comes from the stoichiometry. And 09 and 10 must agree with each
other, because the two solvers receive an identical linear program and share
nothing else.

---

## Why each folder carries its own kinetics headers

CompLaB compiles its rate laws in. `defineKinetics.hh` and
`defineAbioticKinetics.hh` are `#include`d, not read at run time, so two cases
with different chemistry need different binaries. `runAllExamples.sh` copies
each case's headers in and rebuilds, restoring your own on exit — including on
Ctrl-C.

That is not this suite being awkward; it is how the code is built. The 2D suite
worked the same way: every `example/scalability` case shipped its own
`defineKinetics.hh`.

The shipped headers implement the 95-substrate uranium network. Dropped into a
2-substrate case they index `C[0]` to `C[94]` and read past the end of the
array, so every example replaces them. `validateExamples.py` checks each
header's own `C.size()` guard against the case's substrate count, which catches
exactly the mistake of forgetting to.

Three cases need a different cmake configuration: 09 and 12 need
`-DENABLE_GLPK=ON`, 10 needs `-DENABLE_COBRAPY=ON`. The runner groups them, so
the tree is reconfigured three times rather than sixteen.

---

## validateExamples.py

A one-second structural check, run before you build anything. It verifies:

- every XML is well formed, and **no comment contains a double hyphen** — XML
  forbids it and this codebase has tripped over it repeatedly
- every tag name appears in the set the **source actually reads**, extracted
  mechanically from CompLB3D rather than typed out. A tag renamed upstream
  shows up as an unknown tag instead of as silence
- the vector-length and presence rules that terminate a run: `Kc` and `Vmax`
  lengths against the substrate count, `initial_densities` against the material
  numbers, `viscosity_ratio_in_biofilm` present for exactly the seeded microbes,
  `biomass_diffusion_coefficients` present when the solver is FD, FD rejected
  for planktonic microbes
- reaction types spelled acceptably, and the matching `enable_*` switch on
- geometry files: the right voxel count, and every material number the XML
  refers to actually present
- metabolic model files: `nmet`, `nrxn`, and the lengths of S, b, c, lb, ub;
  exchange indices inside range
- precipitation and dissolution pointing at substrates that are actually
  `<immobile>true</immobile>`
- the kinetics headers **compile**, at `-Wall -Wextra`

It is not Palabos. It builds no lattices and runs nothing.

**It was checked against deliberate breakage**, which is the only way to trust
a checker. Ten mutations were introduced one at a time — a short `Kc` vector, an
exchange index past the end of the model, `gplk` for `glpk`, the enable switch
turned off, a mistyped tag, FD on a planktonic microbe, a truncated S matrix, a
double hyphen in a comment, a kinetics guard written for other chemistry, a
syntax error in a rate law. **All ten were caught.**

---

## Extending the suite

`makeExamples.py` generates the folders from shared fragments in `_build/`.
Edit and regenerate rather than hand-patching sixteen XMLs; that is how the
suite stays consistent after the validator finds something.

A case is a `case(...)` call: a number, a name, the XML, a README, a geometry,
and optionally its own kinetics headers and a metabolic model.

---

## Not covered here

- **checkpoint restart** — `<read_NS_file>` and `<read_ADE_file>`, which resume
  from a `.chk` written by `<save_CHK_interval>`. Set `save_CHK_interval` in any
  case above, run it, then rerun with the read flags on.
- **MPI** — every case runs under `mpirun` unchanged; the domain is too small
  for the scaling to mean anything.
- **the scalability benchmarks** from the 2D suite (`example/scalability`,
  four cases). Those are performance measurements rather than capability
  demonstrations, and are still to be ported.

---

## What has actually been run

The suite has been built and run against Palabos v2.3.0. Cases 01 to 12 and 16
pass. Doing it found four real bugs, listed in `WHAT_TO_UPLOAD.md` at the top of
the repository; two of them were in the shipped code and one made every example
in this folder fail to compile.

**Cases 13, 14 and 15 report negative concentrations and are marked FAIL.** The
cause is not the chemistry: with the abiotic kinetics *and* the precipitation
solver both switched off, Fe²⁺ still reaches exactly −0.18889 at iteration 0.
It is an undershoot of the ADE Dirichlet inlet on the advected species, in the
transport core, unaffected by Peclet and by tau. It is recorded here rather than
papered over.

**Case 10 needs COBRApy** in the same interpreter the executable was linked
against, which is not the same thing as the `python3` on your PATH. If it fails
to import `cobra`, the run now prints which interpreter is embedded.
