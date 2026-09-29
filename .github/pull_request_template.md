## Problem

<!-- What problem does this PR solve? -->

## Changes

<!-- Summarize the smallest coherent implementation change. -->

## Architecture / reuse

<!-- What existing modules/helpers were reused or extended? -->

## Validation

- [ ] Targeted tests pass
- [ ] `pytest tests/ -v` passes
- [ ] No duplicate implementation was introduced
- [ ] CLI/docs updated if user-visible behavior changed

## Scientific invariants

- [ ] Tree topology semantics preserved
- [ ] Trait values are not silently altered/imputed
- [ ] Time interpretation is not assumed without justified branch-length semantics
- [ ] No dataset-specific taxa/node IDs/values were hard-coded
