# Example 07: biotic kinetics, LBM biomass solver

One organism, Monod growth on one donor, biomass transported by LBM.
Examples 05, 06 and 07 differ only in `<solver_type>`.

## What to expect

- biomass rising from the seeded patch, then levelling off as
  growth and decay balance
- the donor drawn down around the colony and replenished from the
  left boundary
- a product plume spreading from the colony
- no `[NEG!]` warnings: the donor must never go negative

## The switches this case exercises

`<solver_type>LBM</solver_type>`, lattice Boltzmann: biomass gets its own advection diffusion
     lattice, the same machinery the solutes use

## Notes

LBM costs a full extra lattice per microbe. It is the right choice
when biomass genuinely moves with the flow, and overkill when it
does not.

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
