# Example 06: biotic kinetics, FD biomass solver

One organism, Monod growth on one donor, biomass transported by FD.
Examples 05, 06 and 07 differ only in `<solver_type>`.

## What to expect

- biomass rising from the seeded patch, then levelling off as
  growth and decay balance
- the donor drawn down around the colony and replenished from the
  left boundary
- a product plume spreading from the colony
- no `[NEG!]` warnings: the donor must never go negative

## The switches this case exercises

`<solver_type>FD</solver_type>`, finite difference: biomass diffuses with its own coefficient,
     solved on the same grid

## Notes

FD is the only solver that requires `<biomass_diffusion_coefficients>`,
and it is rejected outright for planktonic microbes: a microbe with
no `<microbeN>` entry under `material_numbers` cannot use it.

## Running it

The rate laws are compiled in, so each case needs its own build.

```bash
cp defineKinetics.hh defineAbioticKinetics.hh  <path to CompLB3D>/
cd <path to CompLB3D> && mkdir -p build && cd build
cmake -DPALABOS_ROOT=<path to palabos-v2.3.0> .. && make -j
cd .. && cp <path to this folder>/CompLaB.xml .
cp -r <path to this folder>/input .
./complab
```

`runAllExamples.sh` in the parent folder does all of that for every case, and
groups the cases by cmake configuration so the tree is reconfigured three times
rather than sixteen.
