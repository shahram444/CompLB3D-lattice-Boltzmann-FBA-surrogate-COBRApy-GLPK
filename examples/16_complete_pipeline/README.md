# Example 16: the whole pipeline, from the XML alone

The four things that used to be manual, all in one file. The metabolic
model comes out of the bundle; the exchange reactions are named rather
than numbered; the surrogate network is fitted during start-up from
that model; and the run writes its own summary CSV.

There is no preparation step. No extractMM.py, no training scripts, no
post-processing to find out whether anything was conserved.

## What to expect

- `[MODEL] unpacking bundled models/e_coli_core.xml.gz`, then
  `revision check: 72 species, 95 reaction tags, as the manifest expects`
- `resolved 3 exchange reaction name(s) against e_coli_core`
- `[SRG] no weights file; training one now. This happens once.`
  followed by a few hundred linear programs, then a held-out R^2 above
  0.99
- **`output/ecoli.srg` written.** Run the case again: it loads instead
  of training, and start-up is instant
- `output/summary.csv`, one row per 100 steps
- at the end, `[SRG] runtime network: N evaluation(s), 0 clamped`. A
  non-zero clamp count means the simulation left the box the network
  was fitted on

## The switches this case exercises

`<model_source>`, `<exchange_reaction_names>`, `<surrogate>` with
`<train_if_missing>`, and `<diagnostics>`. Needs `-DENABLE_GLPK=ON`,
because training sweeps a linear program.

## Notes

**This run has no GLPK microbe**, so there is no persistent linear
program to sweep. A temporary one is built from this microbe's own
model for training and released as soon as the fit is done. The log
says so.

**Why the geometry is still a file here.** `<generate>` in `<domain>`
can build the pore space from the XML, and example 01 shows the tags
commented in place. It writes pore, solid and wall only -- it cannot
seed the microbe material numbers this case needs, so the seeded
geometry is shipped.

**Comparing against the FBA it replaces.** Set `<enable_fba_glpk>true`
and `<reaction_type>glpk</reaction_type>`, delete the `<surrogate>`
block, and rerun. On the reference run the two agreed to four parts in
ten thousand in final biomass, and the surrogate was about 150 times
faster.

## Running it

The rate laws are compiled in, so each case needs its own build.

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && mkdir -p build && cd build
cmake -DPALABOS_ROOT=<path to palabos-v2.3.0> -DENABLE_GLPK=ON .. && make -j
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
# <model_source> looks for the bundle in the WORKING directory,
# and models/ is already at the top of the CompLB3D tree, so running
# from there is all this needs.
./complab
```

**`-DENABLE_GLPK=ON` is not optional for this case.** Without it the run
stops at start-up naming the flag.

`runAllExamples.sh` in the parent folder does all of that for every case, and
groups the cases by cmake configuration so the tree is reconfigured three times
rather than sixteen.
