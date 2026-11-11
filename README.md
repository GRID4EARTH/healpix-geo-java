# healpix-geo-java

HEALPix on the WGS84 ellipsoid in Java: the NESTED cell id of a longitude/latitude, and the
centre of a cell. It follows the ellipsoidal HEALPix definition of the
[CF conventions 1.13](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.13/cf-conventions.html#healpix) (grid mapping `healpix`): geodetic latitude is
converted to authalic latitude, then standard HEALPix is applied.

- Cell ids: NESTED indexing (CF `indexing_scheme: nested`), as a `long`.
- Resolution: `level` 0 to 29 (CF `refinement_level`), 12 × 4^level cells. The healpix-geo
  `nested` functions call the same number `depth`.

It is a small port of two functions of [healpix-geo](https://github.com/GRID4EARTH/healpix-geo),
the library [GRID4EARTH](https://grid4earth.eu/) uses for HEALPix on the ellipsoid, for Java
environments such as GBIF's occurrence cubes. No dependencies.

```java
import eu.grid4earth.healpix.HealpixWgs84;

long cell = HealpixWgs84.cell(10, 2.35, 48.85);       // level, longitude, latitude (degrees)
double[] centre = HealpixWgs84.centre(10, cell);      // {longitude in [0, 360[, latitude}
```

## Agreement with healpix-geo

The tests compare both functions with values written by healpix-geo itself, at every level,
including cell corners, the poles and the polar-cap transition latitude. A weekly CI job repeats
the comparison with the latest healpix-geo release. Details, and what was ported from where, are
in [UPSTREAM.md](UPSTREAM.md).

Cell ids are identical except for points lying on a cell corner, or within about 1e-12 degree
of one (0.1 micrometre), where a last-bit difference in `sin`/`cos` between Java and Rust may
select the neighbouring cell. Cell centres agree within 3e-14 degree.

## Build

```bash
mvn test
```

Requires Java 17 or later.

## Licence

Apache License 2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
