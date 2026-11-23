# Frozen API migration reference

`rsm-api-b0d059c.R` was captured from package commit
`b0d059c0a4a985b9f623a663c2329965ad91e84e` before the argument-based API change.
It contains synthetic draws (seed 37, three replications), Stan input, and
generated Stan code for the historical defaults, a custom raw-scale prior
including Student-t, and the previous explicit marginal Gaussian specification.
The corresponding inputs are in `test-api-compatibility.R`.

Hexadecimal numeric literals preserve double precision exactly; matrix
dimensions and names are preserved. These are regression references, not
estimates of sampling accuracy. Independent density/gradient verification is
performed separately by `tools/check-rsm-stan.R`.
