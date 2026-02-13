# Developer Guide

## Make Commands

This repository defines the following `make` targets:

- `make tests`: run the package test suite (`test/runtests.jl`).
- `make doc`: build the documentation site from `docs/`.
- `make format`: format all tracked `*.jl` files with Runic.

## Consistent Formatting (Runic)

Corpuscles.jl uses [Runic](https://github.com/fredrikekre/Runic.jl) for consistent Julia formatting.

- CI validates formatting on pull requests.
- Locally, run `make format` before committing.

To run a formatting check without modifying files:

```sh
git ls-files -z -- '*.jl' | xargs -0 --no-run-if-empty julia --project=@runic --startup-file=no -e 'using Runic; exit(Runic.main(ARGS))' -- --check --diff
```
