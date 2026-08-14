# Example 12: combined reaction types, FBA plus hand written kinetics

The example 09 organism with `<reaction_type>glpk_and_kinetics`, so
the linear program and `defineKinetics.hh` both run and their rates
add. The kinetics file contributes maintenance: a donor drain that
yields no growth.

## What to expect

- the same growth rate as example 09, near 8 mmol/gDW/h, because
  maintenance does not add biomass
- the donor drawn down **faster** than in example 09, by exactly the
  maintenance term
- the metabolic summary showing 1 GLPK microbe, and the kinetics
  path also active

## The switches this case exercises

`<reaction_type>glpk_and_kinetics</reaction_type>`, both
`<enable_kinetics>` and `<enable_fba_glpk>` true, and both
`<maximum_uptake_flux>` and `<fba_maximum_uptake_flux>` present.

## Notes

Diff this case's donor profile against example 09's. The difference
IS the maintenance term, and it is the cheapest way to confirm the
two paths are genuinely adding rather than one silently winning.

## Running it

The rate laws are compiled in, so each case needs its own build.

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && mkdir -p build && cd build
cmake -DPALABOS_ROOT=<path to palabos-v2.3.0> -DENABLE_GLPK=ON .. && make -j
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

**`-DENABLE_GLPK=ON` is not optional for this case.** Without it the run
stops at start-up naming the flag.

`runAllExamples.sh` in the parent folder does all of that for every case, and
groups the cases by cmake configuration so the tree is reconfigured three times
rather than sixteen.
