## Resubmission of an archived package

This is a resubmission of `bayeslm`, which is currently archived on CRAN.

Since the archive, the package has been updated to address CRAN-facing issues and
to improve the C++ implementation used by the samplers. In particular:

- the package now caps `RcppParallel` execution at no more than 2 threads,
  matching CRAN policy for packages that use parallel computation;
- package metadata and repository links were refreshed;
- the package was rechecked locally with current R and external tools.

## R CMD check results

I ran `R CMD check --as-cran` with:

- R 4.5.3 on macOS
- `qpdf` installed
- HTML Tidy installed
- the `V8` package installed

Result:

- 0 errors
- 0 warnings
- 1 note

The remaining note is the standard CRAN incoming note that this is a new
submission of a package that is currently archived:

```text
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jingyu He <jingyuhe@cityu.edu.hk>'

New submission

Package was archived on CRAN

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2025-11-27 as issues were not corrected
    despite reminders.
```

## Additional checks

I also submitted the package tarball to Win-builder for `R-devel` on Windows.
