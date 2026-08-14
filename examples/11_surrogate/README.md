# Example 11: surrogate model, a trained network instead of FBA

The shipped Geobacter network standing in for a flux balance solve.
Two inputs, acetate and Fe(III) uptake; one output, growth rate.

## What to expect

- growth rates in the low 1e-2 1/h range: the shipped network's whole
  output range is 0 to 0.0527 1/h
- a step time far below examples 09 and 10 on the same domain, which
  is the entire reason to use a surrogate
- both substrates drawn down together

## The switches this case exercises

`<enable_surrogate>true`, `<reaction_type>surrogate</reaction_type>`,
and `surrogateModel.hh` compiled in.

## Notes

**The Vmax values matter here in a way they do not elsewhere.** The
shipped network has never seen a Fe(III) uptake above 0.5 mmol/gDW/h.
`<fba_maximum_uptake_flux>` is set to 0.4 to stay inside that. Raise
it and the run leaves the training box, where a network extrapolates
badly and returns confident nonsense.

**There are now two ways to get a network in, and this case shows the
older one.** Here the weights are compiled into `surrogateModel.hh`, so
changing them means editing that file and rebuilding.

The alternative is `<surrogate><weights_file>`, which reads the weights
at run time and can also FIT them during start-up from a metabolic
model. **Example 16 shows that.** Both routes are supported, and a run
with no `<surrogate>` block behaves exactly as this one does.

To fit your own network offline, see `surrogate_training/` in the parent
folder. `inspectSurrogate.py` prints any network's valid range,
recovered by inverting the mapminmax scaling.

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
