"""Reference values from healpix-geo, which the Java tests must reproduce.

    python reference/make_reference.py [output directory]

Writes, gzipped, to src/test/resources/reference/ by default:

- cells.csv.gz:   kind, depth, lon, lat, cell   (healpix-geo lonlat_to_healpix, WGS84)
- centres.csv.gz: depth, cell, lon, lat         (healpix-geo healpix_to_lonlat, WGS84)
- healpix-geo-version.txt: the healpix-geo version that produced them

Floats are written with repr(), so they read back exactly.
"""

import gzip
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
from healpix_geo import nested

ELLIPSOID = "WGS84"
DEPTHS = range(30)
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parents[1] / "src/test/resources/reference"
rng = np.random.default_rng(20261104)


def cells_rows():
    rows = []

    def add(kind, depth, lon, lat):
        lon = np.asarray(lon, dtype=float).ravel()
        lat = np.clip(np.asarray(lat, dtype=float).ravel(), -90.0, 90.0)
        ids = nested.lonlat_to_healpix(lon, lat, depth, ellipsoid=ELLIPSOID)
        rows.extend(f"{kind},{depth},{float(a)!r},{float(b)!r},{int(c)}" for a, b, c in zip(lon, lat, ids))

    for d in DEPTHS:
        # uniform on the sphere; longitudes in [-180, 540) exercise the wrap
        add("random", d, rng.uniform(-180, 540, 4000), np.degrees(np.arcsin(rng.uniform(-1, 1, 4000))))
        cells = rng.integers(0, 12 * 4**d, 300, dtype=np.uint64)
        clon, clat = nested.healpix_to_lonlat(cells, d, ellipsoid=ELLIPSOID)
        add("centre", d, clon, clat)
        vlon, vlat = nested.vertices(cells, d, ellipsoid=ELLIPSOID)
        vlon, vlat = np.asarray(vlon).ravel(), np.asarray(vlat).ravel()
        add("vertex", d, vlon, vlat)
        for eps in (1e-12, 1e-9, 1e-6):
            jitter = rng.choice([-eps, eps], size=(2, vlon.size))
            add(f"near{eps:g}", d, vlon + jitter[0], vlat + jitter[1])

    t = 41.937853910160136  # geodetic latitude of the polar-cap transition (asin(2/3) authalic)
    lats = [90.0, -90.0, 0.0, -0.0, t, -t, np.nextafter(t, 0), np.nextafter(t, 90), 89.999999999, -89.999999999]
    lons = [0.0, -0.0, 360.0, -360.0, 720.0, 45.0, 90.0, 135.0, 180.0, -180.0, 359.99999999999994, -1e-300, 1e-300]
    lo, la = np.meshgrid(lons, lats)
    for d in DEPTHS:
        add("special", d, lo, la)
    return "kind,depth,lon,lat,cell", rows


def centres_rows():
    rows = []
    for d in DEPTHS:
        n = 12 * 4**d
        cells = np.unique(np.concatenate([
            rng.integers(0, n, min(3000, n), dtype=np.uint64),
            np.array([0, n - 1, n // 2, 4 * 4**d, 8 * 4**d - 1], dtype=np.uint64),
        ]))
        lon, lat = nested.healpix_to_lonlat(cells, d, ellipsoid=ELLIPSOID)
        rows.extend(f"{d},{int(c)},{float(a)!r},{float(b)!r}" for c, a, b in zip(cells, lon, lat))
    return "depth,cell,lon,lat", rows


OUT.mkdir(parents=True, exist_ok=True)
for name, (header, rows) in (("cells.csv.gz", cells_rows()), ("centres.csv.gz", centres_rows())):
    with gzip.open(OUT / name, "wt", compresslevel=9) as f:
        f.write(header + "\n" + "\n".join(rows) + "\n")
    print(f"{name}: {len(rows)} rows")
(OUT / "healpix-geo-version.txt").write_text(version("healpix-geo") + "\n")
print("healpix-geo", version("healpix-geo"))
