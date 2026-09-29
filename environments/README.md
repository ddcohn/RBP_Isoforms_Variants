# Conda environments

Exact, exported package specifications for the four conda environments
this project's pipeline depends on. Before this, these environments
existed only as already-built directories on Hoffman2's filesystem --
nothing about what packages/versions they actually contained was
recorded anywhere, so nobody besides the person who originally built
them could recreate one from scratch. These `.yml` files fix that.

Exported with `conda env export --no-builds` (drops OS-specific build
hashes, which are the most common cause of an exported environment
failing to install on a different machine or even a different Hoffman2
node type; the exact version pins are kept, so this is still a precise
reproduction, not a loose one) and with the machine-specific `prefix:`
line stripped (an absolute path on one person's home directory, not
needed to recreate the environment elsewhere).

| File | Used for | Where |
|---|---|---|
| `spliceai.yml` | SpliceAI splice-effect prediction | `12_variant_effect_prediction/splicing/` |
| `deeploc2.yml` | DeepLoc 2.1 subcellular localization prediction | `12_variant_effect_prediction/localization_condensate/` |
| `protgps.yml` | protGPS condensate-localization prediction | `12_variant_effect_prediction/localization_condensate/` |
| `snakemake.yml` | Orchestrates the whole pipeline (`workflow/Snakefile`) | `12_variant_effect_prediction/workflow/` |

## Recreating an environment

```
module load miniforge/23.11.0
conda env create -f environments/spliceai.yml
```

(substitute the other three the same way). This is also exactly how to
set up the pipeline on a fresh Hoffman2 account, or in principle any
other Linux machine with conda/mamba installed -- see the top-level
`12_variant_effect_prediction/README.md`'s "Reproduce from scratch"
section for the rest of the setup.

## Keeping these current

These are a snapshot as of when they were exported (2026-09-29) --
if a package gets upgraded in one of these environments later, these
files won't update themselves. Re-export with the same command above
and commit the new version if that happens.
