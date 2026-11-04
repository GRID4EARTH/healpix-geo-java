package eu.grid4earth.healpix;

/**
 * HEALPix on the WGS84 ellipsoid, NESTED cell ids, as in OGC API - DGGS Annex B.8
 * (https://www.opengis.net/def/dggrs/OGC/1.0/HEALPix): geodetic latitude converted to
 * authalic latitude, then standard HEALPix.
 *
 * <p>A port of two functions of healpix-geo (https://github.com/GRID4EARTH/healpix-geo,
 * Apache-2.0), kept in agreement with it by the tests and by a weekly CI job (see UPSTREAM.md):
 *
 * <ul>
 *   <li>{@link #cell}: healpix-geo {@code nested::coordinates::lonlat_to_healpix}
 *   <li>{@link #centre}: healpix-geo {@code nested::coordinates::healpix_to_lonlat}
 * </ul>
 *
 * <p>which rely on:
 *
 * <ul>
 *   <li>geodetic/authalic latitude: geodesy 0.15 {@code ellipsoid/latitudes.rs}, coefficients
 *       {@code AUTHALIC} in {@code ellipsoid/constants.rs} (Karney 2022, eq. A19 and A20);
 *   <li>NESTED hash and cell centre: cdshealpix 0.9.1 {@code nested::Layer::hash},
 *       {@code nested::Layer::center} and {@code unproj} (CDS, Apache-2.0 OR MIT).
 * </ul>
 */
public final class HealpixWgs84 {

  // WGS84 (EPSG:7030)
  private static final double FLATTENING = 1.0 / 298.257223563;

  // Karney (2022) eq. A19: geodetic -> authalic, polynomials in the third flattening n
  private static final double[][] AUTHALIC_FWD = {
    {-4. / 3., -4. / 45., 88. / 315., 538. / 4725., 20824. / 467775., -44732. / 2837835.},
    {0., 34. / 45., 8. / 105., -2482. / 14175., -37192. / 467775., -12467764. / 212837625.},
    {0., 0., -1532. / 2835., -898. / 14175., 54968. / 467775., 100320856. / 1915538625.},
    {0., 0., 0., 6007. / 14175., 24496. / 467775., -5884124. / 70945875.},
    {0., 0., 0., 0., -23356. / 66825., -839792. / 19348875.},
    {0., 0., 0., 0., 0., 570284222. / 1915538625.},
  };

  // Karney (2022) eq. A20: authalic -> geodetic
  private static final double[][] AUTHALIC_INV = {
    {4. / 3., 4. / 45., -16. / 35., -2582. / 14175., 60136. / 467775., 28112932. / 212837625.},
    {0., 46. / 45., 152. / 945., -11966. / 14175., -21016. / 51975., 251310128. / 638512875.},
    {0., 0., 3044. / 2835., 3802. / 14175., -94388. / 66825., -8797648. / 10945935.},
    {0., 0., 0., 6059. / 4725., 41072. / 93555., -1472637812. / 638512875.},
    {0., 0., 0., 0., 768272. / 467775., 455935736. / 638512875.},
    {0., 0., 0., 0., 0., 4210684958. / 1915538625.},
  };

  private static final double[] FOURIER_FWD = fourierCoefficients(AUTHALIC_FWD);
  private static final double[] FOURIER_INV = fourierCoefficients(AUTHALIC_INV);

  // cdshealpix constants
  private static final double TRANSITION_LATITUDE = 0.729_727_656_226_966_3; // asin(2/3)
  private static final double ONE_OVER_TRANSITION_Z = 1.5;
  private static final double SQRT6 = 2.449_489_742_783_178;
  private static final double FOUR_OVER_PI = 4.0 / Math.PI;
  private static final double FRAC_PI_4 = Math.PI / 4.0;
  private static final double HALF = 0.5;
  private static final long SIGN_BIT = 0x8000000000000000L;

  private HealpixWgs84() {}

  /** NESTED cell id at {@code depth} (0..29) of a WGS84 point given in degrees. */
  public static long cell(int depth, double lonDeg, double latDeg) {
    if (depth < 0 || depth > 29) {
      throw new IllegalArgumentException("depth must be between 0 and 29");
    }
    if (!(latDeg >= -90.0 && latDeg <= 90.0)) {
      throw new IllegalArgumentException("latitude must be in [-90, 90]");
    }
    double lon = Math.toRadians(remEuclid(lonDeg, 360.0));
    double lat = geographicToAuthalic(Math.toRadians(latDeg));
    return hash(depth, lon, lat);
  }

  // ---- authalic latitude (geodesy crate) ----

  static double geographicToAuthalic(double lat) {
    return lat + fourierSin(2.0 * lat, FOURIER_FWD);
  }

  static double authalicToGeographic(double lat) {
    return lat + fourierSin(2.0 * lat, FOURIER_INV);
  }

  private static double[] fourierCoefficients(double[][] poly) {
    double n = FLATTENING / (2.0 - FLATTENING); // third flattening
    double[] out = new double[6];
    for (int i = 0; i < 6; i++) {
      out[i] = n * horner(n, poly[i]);
    }
    return out;
  }

  /**
   * Centre of a NESTED cell at {@code depth}, as {longitude in [0, 360[, WGS84 geodetic latitude},
   * in degrees.
   */
  public static double[] centre(int depth, long cell) {
    if (depth < 0 || depth > 29) {
      throw new IllegalArgumentException("depth must be between 0 and 29");
    }
    if (cell < 0 || cell >= 12L << (2 * depth)) {
      throw new IllegalArgumentException("cell id out of range for this depth");
    }
    double[] lonlat = centreRadians(depth, cell);
    return new double[] {
      remEuclid(Math.toDegrees(lonlat[0]), 360.0),
      Math.toDegrees(authalicToGeographic(lonlat[1]))
    };
  }

  /** Sum c_i x^i by Horner's scheme, with fused multiply-add as Rust's mul_add. */
  private static double horner(double x, double[] c) {
    double v = c[c.length - 1];
    for (int k = c.length - 2; k >= 0; k--) {
      v = Math.fma(v, x, c[k]);
    }
    return v;
  }

  /** Sum c_i sin(i * arg), i = 1..6, by Clenshaw summation. */
  private static double fourierSin(double arg, double[] c) {
    double x = 2.0 * Math.cos(arg);
    double c0 = 0.0, c1 = 0.0;
    for (int k = c.length - 1; k >= 0; k--) {
      double next = Math.fma(x, c0, c[k] - c1);
      c1 = c0;
      c0 = next;
    }
    return Math.sin(arg) * c0;
  }

  private static double remEuclid(double x, double m) {
    double r = x % m;
    return r < 0.0 ? r + Math.abs(m) : r;
  }

  // ---- NESTED hash (cdshealpix Layer::hash) ----

  static long hash(int depth, double lon, double lat) {
    long nside = 1L << depth;
    long timeHalfNside = ((long) depth - 1) << 52;

    double xPm1;
    int q;
    long lonBits = Double.doubleToRawLongBits(lon);
    double lonAbs = Double.longBitsToDouble(lonBits & ~SIGN_BIT);
    double x = lonAbs * FOUR_OVER_PI;
    int qq = saturatingU8(x) | 1;
    if ((lonBits & SIGN_BIT) == 0) {
      xPm1 = x - qq;
      q = (qq & 7) >> 1;
    } else {
      xPm1 = qq - x;
      q = 3 - ((qq & 7) >> 1);
    }

    int d0h;
    double lInD0c, hInD0c;
    if (lat > TRANSITION_LATITUDE) {
      double s = SQRT6 * Math.cos(HALF * lat + FRAC_PI_4);
      d0h = q;
      lInD0c = xPm1 * s;
      hInD0c = 2.0 - s;
    } else if (lat < -TRANSITION_LATITUDE) {
      double s = SQRT6 * Math.cos(HALF * lat - FRAC_PI_4);
      d0h = q + 8;
      lInD0c = xPm1 * s;
      hInD0c = s;
    } else {
      double yPm1 = Math.sin(lat) * ONE_OVER_TRANSITION_Z;
      int q01 = xPm1 > yPm1 ? 1 : 0;
      int q12 = xPm1 >= -yPm1 ? 1 : 0;
      int q1 = q01 & q12;
      int q013 = q01 + (1 - q12);
      lInD0c = xPm1 - ((q01 + q12) - 1);
      hInD0c = yPm1 + q013;
      d0h = (q013 << 2) + ((q + q1) & 3);
    }

    long i = scaledU32(timeHalfNside, hInD0c + lInD0c);
    long j = scaledU32(timeHalfNside, hInD0c - lInD0c);
    if (i == nside) i = nside - 1;
    if (j == nside) j = nside - 1;
    return ((long) d0h << (2 * depth)) | interleave(i) | (interleave(j) << 1);
  }

  /** f64::from_bits(timeHalfNside + v.to_bits()) as u32: v * nside/2, saturating cast. */
  private static long scaledU32(long timeHalfNside, double v) {
    double d = Double.longBitsToDouble(timeHalfNside + Double.doubleToRawLongBits(v));
    if (!(d > 0.0)) return 0L; // NaN and negatives saturate to 0
    if (d >= 4294967295.0) return 4294967295L;
    return (long) d;
  }

  private static int saturatingU8(double x) {
    if (!(x > 0.0)) return 0;
    if (x >= 255.0) return 255;
    return (int) x;
  }

  // ---- cell centre (cdshealpix Layer::center + unproj) ----

  private static final double ONE_OVER_SQRT6 = 0.408_248_290_463_863;
  private static final double TRANSITION_Z = 2.0 / 3.0;
  private static final double EPS_POLE = 1e-13;

  static double[] centreRadians(int depth, long hash) {
    int twiceDepth = 2 * depth;
    long nsideMinus1 = (1L << depth) - 1;
    double oneOverNside = Double.longBitsToDouble(Double.doubleToRawLongBits(1.0) - ((long) depth << 52));
    int d0h = (int) (hash >>> twiceDepth);
    long ij = hash & ((1L << twiceDepth) - 1);
    int i = (int) compact(ij);
    int j = (int) compact(ij >>> 1);
    // rotate45_scale2, then shift from small-cell centre to base-cell centre
    int h = i - j;
    int l = (i + j) - (int) nsideMinus1;
    double x = h * oneOverNside;
    double y = l * oneOverNside;
    int offsetY = 1 - (d0h >> 2);
    int offsetX = ((d0h & 3) << 1) | (offsetY & 1);
    x += offsetX;
    y += offsetY;
    x += (double) ((Double.doubleToRawLongBits(x) & SIGN_BIT) >>> 60); // x < 0 -> x + 8
    return unproj(x, y);
  }

  private static double[] unproj(double x, double y) {
    long xBits = Double.doubleToRawLongBits(x), yBits = Double.doubleToRawLongBits(y);
    long xSign = xBits & SIGN_BIT, ySign = yBits & SIGN_BIT;
    double xAbs = Double.longBitsToDouble(xBits & ~SIGN_BIT);
    double yAbs = Double.longBitsToDouble(yBits & ~SIGN_BIT);
    int oddFloor = saturatingU8(xAbs) | 1;
    int offset = oddFloor & 7;
    double lon = xAbs - oddFloor;
    double lat = yAbs;
    if (lat <= 1.0) {
      lat = Math.asin(lat * TRANSITION_Z);
    } else {
      lat = 2.0 - lat;
      if (lat > EPS_POLE) {
        lon /= lat;
        lon = Math.max(-1.0, Math.min(1.0, lon));
      }
      lat *= ONE_OVER_SQRT6;
      lat = 2.0 * Math.acos(lat) - Math.PI / 2.0;
    }
    lon += offset;
    lon = Double.longBitsToDouble(Double.doubleToRawLongBits(lon) | xSign);
    lat = Double.longBitsToDouble(Double.doubleToRawLongBits(lat) | ySign);
    return new double[] {lon * FRAC_PI_4, lat};
  }

  /** Gather the even bits of v into the low 32 bits (inverse of interleave). */
  private static long compact(long v) {
    v &= 0x5555555555555555L;
    v = (v | (v >>> 1)) & 0x3333333333333333L;
    v = (v | (v >>> 2)) & 0x0F0F0F0F0F0F0F0FL;
    v = (v | (v >>> 4)) & 0x00FF00FF00FF00FFL;
    v = (v | (v >>> 8)) & 0x0000FFFF0000FFFFL;
    v = (v | (v >>> 16)) & 0x00000000FFFFFFFFL;
    return v;
  }

  /** Spread the low 32 bits of v to the even bit positions. */
  private static long interleave(long v) {
    v &= 0xFFFFFFFFL;
    v = (v | (v << 16)) & 0x0000FFFF0000FFFFL;
    v = (v | (v << 8)) & 0x00FF00FF00FF00FFL;
    v = (v | (v << 4)) & 0x0F0F0F0F0F0F0F0FL;
    v = (v | (v << 2)) & 0x3333333333333333L;
    v = (v | (v << 1)) & 0x5555555555555555L;
    return v;
  }
}
