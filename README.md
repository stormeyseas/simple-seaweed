# simple-seaweed

A Shiny app that simulates seaweed growth and nitrogen removal for several
species under user-chosen environmental conditions, using the
[macrogrow](https://github.com/stormeyseas/macrogrow) model.

## Running locally

Open the project (renv will restore packages) and run `shiny::runApp()`.

## Deployment

The app is published as a [Shinylive](https://posit-dev.github.io/r-shinylive/)
site on GitHub Pages: R runs in the visitor's browser via webR, so no server
is needed. Deployment is automated by
[`.github/workflows/deploy-shinylive.yml`](.github/workflows/deploy-shinylive.yml)
on every push to `main`.

Because `macrogrow` is not on CRAN, webR needs a pre-built WebAssembly binary
of it. The macrogrow repo builds one automatically whenever a GitHub release
is published there (see its `release-wasm.yml` workflow), and this repo's
workflow installs `macrogrow` at the release tag set in `MACROGROW_REF`.

To update macrogrow in the deployed app: publish a new release in the
macrogrow repo, then bump `MACROGROW_REF` in the workflow here.

Note: the app avoids the `units` package (plain conversion factors live in
`R/00_conversions.R`) to keep the browser bundle small.
