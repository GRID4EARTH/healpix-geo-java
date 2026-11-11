# Keeping in agreement with healpix-geo

This library is a port, so it can drift when healpix-geo or the Rust crates beneath it change.
Agreement is checked, not assumed.

| Ported here | Upstream | Version ported |
|---|---|---|
| `HealpixWgs84.cell` | healpix-geo `nested::coordinates::lonlat_to_healpix` | healpix-geo 0.4.1 |
| `HealpixWgs84.centre` | healpix-geo `nested::coordinates::healpix_to_lonlat` | healpix-geo 0.4.1 |
| NESTED hash, cell centre | cdshealpix `nested::Layer::hash`, `nested::Layer::center`, `unproj` | cdshealpix 0.9.1 (rev d9d017c, as pinned by healpix-geo 0.4.1) |
| geodetic/authalic latitude | geodesy `ellipsoid/latitudes.rs`, `AUTHALIC` in `ellipsoid/constants.rs` | geodesy 0.15.0 |

## How agreement is checked

1. `reference/make_reference.py` writes reference values with healpix-geo itself:
   random points, cell centres, cell corners exact and nudged by 1e-12, 1e-9 and 1e-6 degree,
   poles, the polar-cap transition latitude and awkward longitudes, at every level 0 to 29.
   The healpix-geo version that produced them is in
   `src/test/resources/reference/healpix-geo-version.txt`.
2. The tests (`mvn test`) require identical cell ids, except points on a cell corner or within
   1e-12 degree of one (about 0.1 micrometre), where a last-bit difference in `sin`/`cos`
   between Java and Rust may pick the neighbouring cell (allowed rate 1e-3; measured 2e-4).
   Cell centres must agree within 1e-12 degree and map back to their own cell.
3. A weekly CI job (`.github/workflows/upstream.yml`) regenerates the reference with the
   latest healpix-geo release and runs the tests against it.

## When healpix-geo changes

If the weekly job fails, compare the upstream functions listed above with the version ported,
port the change, regenerate the committed reference
(`python reference/make_reference.py`), and update the table.
