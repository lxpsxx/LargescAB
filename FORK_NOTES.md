# LargescAB Fork Notes

This package is a fork and large-scale engineering extension of the original
`scAB` package.

## Before publishing

1. Keep the original license terms (`GPL-3`) and attribution.
2. Ensure the package name (`LargescAB`) is consistent across:
   `DESCRIPTION`, `README`, and repository title.
3. Keep public claims precise:
   preserve core scAB methodology, but explicitly state engineering changes for
   scalability.
4. Do not include local run outputs, checkpoints, logs, or private raw data in
   this package repository.
5. Confirm that bundled data files are redistributable.
6. Make sure examples in README are runnable or clearly marked as
   non-executable for large inputs.

## Recommended next checks

1. Run `R CMD build` and `R CMD check` on this clean copy.
2. Confirm the install path:
   `devtools::install_github("<user>/LargescAB")`.
3. Add a short `NEWS.md` entry for your first public release.

