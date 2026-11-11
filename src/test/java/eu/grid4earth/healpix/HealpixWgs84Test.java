package eu.grid4earth.healpix;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import org.junit.jupiter.api.Test;

/**
 * Agreement with healpix-geo, on the reference files written by reference/make_reference.py.
 * The system property {@code reference.dir} points the tests at freshly generated files (the
 * weekly upstream check); otherwise the committed ones are used.
 */
class HealpixWgs84Test {

  /**
   * Points on a cell corner, or within 1e-12 degree of one, may land in a neighbouring cell:
   * a last-bit difference in sin/cos between Java and Rust decides the side. Allowed rate.
   */
  private static final double CORNER_TOLERANCE = 1e-3;

  /** Largest allowed centre difference, in degrees (about 0.1 micrometre). */
  private static final double CENTRE_TOLERANCE_DEG = 1e-12;

  @Test
  void cellIdsMatchHealpixGeo() throws IOException {
    long[] cornerCount = new long[2]; // points, mismatches
    List<String> failures = new ArrayList<>();
    for (String[] p : rows("cells.csv.gz")) {
      String kind = p[0];
      int level = Integer.parseInt(p[1]);
      long expected = Long.parseUnsignedLong(p[4]);
      long got = HealpixWgs84.cell(level, Double.parseDouble(p[2]), Double.parseDouble(p[3]));
      boolean corner = kind.equals("vertex") || kind.equals("near1e-12");
      if (corner) {
        cornerCount[0]++;
        if (got != expected) cornerCount[1]++;
      } else if (got != expected && failures.size() < 10) {
        failures.add(String.join(",", p) + " -> " + got);
      }
    }
    assertTrue(failures.isEmpty(), "cell ids differ from healpix-geo:\n" + String.join("\n", failures));
    double rate = (double) cornerCount[1] / cornerCount[0];
    assertTrue(rate <= CORNER_TOLERANCE,
        cornerCount[1] + " of " + cornerCount[0] + " points on cell corners differ (rate " + rate + ")");
  }

  @Test
  void centresMatchHealpixGeoAndMapBackToTheirCell() throws IOException {
    for (String[] p : rows("centres.csv.gz")) {
      int level = Integer.parseInt(p[0]);
      long cell = Long.parseUnsignedLong(p[1]);
      double[] got = HealpixWgs84.centre(level, cell);
      double dLon = Math.abs(got[0] - Double.parseDouble(p[2]));
      dLon = Math.min(dLon, 360.0 - dLon);
      double dLat = Math.abs(got[1] - Double.parseDouble(p[3]));
      String where = "level " + level + " cell " + cell;
      assertTrue(dLon <= CENTRE_TOLERANCE_DEG && dLat <= CENTRE_TOLERANCE_DEG,
          where + ": centre differs by " + dLon + ", " + dLat + " degrees");
      assertEquals(cell, HealpixWgs84.cell(level, got[0], got[1]), where + ": centre not in its cell");
    }
  }

  @Test
  void refusesInvalidInput() {
    assertThrows(IllegalArgumentException.class, () -> HealpixWgs84.cell(30, 0, 0));
    assertThrows(IllegalArgumentException.class, () -> HealpixWgs84.cell(-1, 0, 0));
    assertThrows(IllegalArgumentException.class, () -> HealpixWgs84.cell(5, 0, 90.5));
    assertThrows(IllegalArgumentException.class, () -> HealpixWgs84.cell(5, 0, Double.NaN));
    assertThrows(IllegalArgumentException.class, () -> HealpixWgs84.centre(0, 12));
    assertThrows(IllegalArgumentException.class, () -> HealpixWgs84.centre(3, -1));
  }

  private static List<String[]> rows(String name) throws IOException {
    String dir = System.getProperty("reference.dir");
    InputStream raw = dir != null
        ? Files.newInputStream(Path.of(dir, name))
        : HealpixWgs84Test.class.getResourceAsStream("/reference/" + name);
    if (raw == null) throw new IOException("missing reference file " + name);
    List<String[]> out = new ArrayList<>();
    try (BufferedReader r = new BufferedReader(
        new InputStreamReader(new GZIPInputStream(raw), StandardCharsets.UTF_8))) {
      r.readLine(); // header
      for (String line; (line = r.readLine()) != null; ) out.add(line.split(","));
    }
    return out;
  }
}
