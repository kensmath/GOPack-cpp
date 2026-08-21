package JNI;

/**
 * JNI bridge to the GOPack C++ core engine (gopack_jni native library).
 *
 * <p>Loads {@code libgopack_jni.dylib} on macOS or {@code gopack_jni.dll} on
 * Windows from the standard java.library.path. Package the appropriate
 * platform binary alongside your application, or call
 * {@link System#load(String)} with an absolute path resolved at runtime if
 * you need to bundle multiple platform binaries in one jar (recommended:
 * ship {@code natives/win64/gopack_jni.dll} and
 * {@code natives/macos-<arch>/libgopack_jni.dylib}, and pick the right one
 * based on {@code os.name} / {@code os.arch} before calling loadLibrary).
 *
 * <p>Every native method below is a direct 1:1 mapping onto the C++ core API
 * in {@code core/include/gopack/Packer.h}; see that header for the
 * authoritative semantics. This class intentionally carries no packing logic
 * of its own -- it is a pure bridge.
 *
 * <p>Package is {@code JNI} (not {@code org.kensmath.gopack}, its name before
 * 8/2026) to match its currently sole consumer's own package layout --
 * CirclePack vendors this file (and {@link GOPackException}/{@link
 * RandomComplexResult}) as plain source under that same package name so they
 * sit alongside its other native bridges. Since JNI resolves native methods
 * by the calling class's fully-qualified name baked into the compiled
 * library's exported symbols, this package doubles as this repository's own
 * package for these files too -- see {@code jni/cpp/JNI_GOPackNative.cpp}'s
 * exported {@code Java_JNI_GOPackNative_*} symbols, renamed to match.
 */
public final class GOPackNative {

    static {
        System.loadLibrary("gopack_jni");
    }

    private GOPackNative() {}

    /**
     * Loads a *.p packing/triangulation file (the CirclePack-compatible
     * format documented in GOPack's docs/GO_Formats.txt) and computes a
     * maximal packing (mode 1: max pack in the disc/plane, or on the sphere
     * if the complex has no boundary).
     *
     * <p>This bridge method only exposes maximal-packing mode (mode 1).
     * Polygonal/rectangle packing (mode 2) is ported in the C++ core
     * ({@code Packer::setMode}/{@code setPolyCenters}/{@code setRectCenters})
     * but not yet wired up to a JNI entry point -- the CLI's {@code
     * --polygon} flag is the only way to reach it today.
     *
     * @param inputPath  path to a *.p triangulation/packing file
     * @param geometryHint reserved for future use (currently ignored --
     *                     geometry is read from the file's GEOMETRY: field);
     *                     pass 0
     * @param tolerance    reserved for future use (currently ignored -- the
     *                     port uses GOPack's fixed 0.01 visual-error cutoff);
     *                     pass 0.0
     * @param maxPasses    upper bound on riffle passes (pass &lt;= 0 to use the
     *                      native default of 200; see continueRiffle's own
     *                      early-exit-on-convergence behavior -- a higher cap
     *                      costs nothing for inputs that converge sooner)
     * @return euclidean radii for every vertex, as an array of length
     *         nodeCount+1 where index v holds vertex v's radius for
     *         v = 1..nodeCount (matching GOPack/CirclePack's own 1-indexed
     *         vertex numbering -- the same convention used throughout the
     *         C++ core); index 0 is unused. Throws {@link GOPackException}
     *         if the native call fails or the file cannot be read.
     */
    public static native double[] computeMaximalPacking(
            String inputPath, int geometryHint, double tolerance, int maxPasses)
            throws GOPackException;

    /**
     * The in-memory counterpart to {@link #computeMaximalPacking}: computes a
     * maximal packing (mode 1) directly from a triangulation already held in
     * memory -- e.g. CirclePack's own per-vertex neighbor ("flower") data --
     * instead of a *.p file path. Skips both the file write a caller would
     * otherwise need to serialize its in-memory complex to text, and the
     * native side's own text parsing (readpack()) to read it back; for large
     * complexes that parsing step can dominate the actual packing
     * computation, which is the whole reason this entry point exists. See
     * {@code gopack::Packer::loadComplex} in {@code core/src/PackerIO.cpp}
     * for the full semantics this mirrors.
     *
     * @param nodeCount number of vertices in the complex
     * @param flowers   length {@code nodeCount+1}, 1-indexed ({@code
     *                  flowers[0]} is ignored/unused). {@code flowers[v]} is
     *                  vertex v's petal list, in the same convention as the
     *                  *.p FLOWERS format: CLOSED (first element == last
     *                  element) iff v is an interior vertex, OPEN (first !=
     *                  last) iff v is a boundary vertex.
     * @param geometry  0 = Euclidean, -1 = Hyperbolic, +1 = Spherical
     *                  (matches the *.p file's GEOMETRY: field and the C++
     *                  core's {@code gopack::Geometry} enum values directly)
     * @param tolerance reserved for future use (currently ignored -- the
     *                  port uses GOPack's fixed 0.01 visual-error cutoff);
     *                  pass 0.0
     * @param maxPasses upper bound on riffle passes (pass &lt;= 0 to use the
     *                  native default of 200)
     * @return euclidean radii for every vertex, in the same {@code
     *         nodeCount+1}-length, 1-indexed convention as {@link
     *         #computeMaximalPacking}. Throws {@link GOPackException} if the
     *         native call fails (e.g. a malformed complex -- wrong array
     *         lengths, no interior vertex, etc.).
     */
    public static native double[] computeMaximalPackingFromComplex(
            int nodeCount, int[][] flowers, int geometry, double tolerance, int maxPasses)
            throws GOPackException;

    // -------------------------------------------------------------------
    // "Raw triangulation" generators -- computeRandomTri/computeRandomSphere/
    // computeRandomRectangle/computeRandomSquare below all back CirclePack's
    // random_tri command, which by design only ever hands back a
    // triangulation (combinatorics + a raw, unpacked layout), never a
    // packing -- whether the requested shape is a sphere, a
    // rectangle/square, or an arbitrary curve-bounded region. None of them
    // call setMode()/riffle() on the native side.
    //
    // Deliberately no "packed"/"riffled" flag anywhere in
    // RandomComplexResult: whether a triangulation later gets packed is a
    // separate, independent action a caller may or may not take afterward
    // (e.g. via computeMaximalPackingFromComplex/
    // computePolygonalPackingFromComplex on the loaded complex) -- it isn't
    // a fact about the triangulation itself, and nothing on either side of
    // the JNI boundary keeps the native Packer around between calls to track
    // it. Each method below still documents what its geometry/corners fields
    // mean (harmless, useful metadata for that later, independent packing
    // step), but none of it should be read as "this has already been
    // packed".
    //
    // See HANDOFFrandomtrinorepack.md for the original request this
    // implements, and computeRandomDisc further below for GOPack-cpp's one
    // deliberately-still-packed random generator (backing the separate
    // random_disc command, which -- unlike random_tri -- genuinely wants a
    // packing).
    // -------------------------------------------------------------------

    /**
     * Generates a random Delaunay triangulation of an arbitrary closed
     * polygonal region -- not just a fixed disc/square/rectangle -- and
     * returns its raw combinatorics and point layout, <b>not</b> a packing.
     * This is the JNI bridge to
     * {@code gopack::Packer::randomTri(intN, bdryN, graph, cent)} with
     * {@code setMode()}/{@code riffle()} deliberately skipped (see {@code
     * jni/cpp/JNI_GOPackNative.cpp}'s
     * {@code Java_JNI_GOPackNative_computeRandomTri} for the full rationale,
     * and the generator-group comment just above this method).
     *
     * <p>This is deliberately not riffled to a packing: for a disc-topology
     * complex, the maximal packing is, by construction/uniformization, the
     * canonical packing that fills the unit disc -- independent of {@code
     * graphXY}'s actual shape, and so would not resemble the input region.
     * Matches CirclePack's own pure-Java {@code RandomTriangulation}/{@code
     * Triangulation} path, which returns an unpacked layout for the same
     * reason.
     *
     * <p>In the returned {@link RandomComplexResult}: {@link
     * RandomComplexResult#radii} is a uniform placeholder ({@code 0.5} for
     * every vertex) with no packing meaning -- nothing has been solved for.
     * {@link RandomComplexResult#centersRe}/{@link
     * RandomComplexResult#centersIm} hold the real output: the actual
     * randomly-placed/Delaunay point positions. {@link
     * RandomComplexResult#geometry} is always Euclidean ({@code 0}), and
     * {@link RandomComplexResult#corners} is always empty (no polygonal mode
     * was entered).
     *
     * @param intN    number of interior points to generate
     * @param bdryN   number of boundary points to generate
     * @param graphXY the closed boundary polygon's vertices, as a flat
     *                x0,y0,x1,y1,... coordinate list (do not repeat the
     *                first point at the end); length must be even and at
     *                least 6 (i.e. at least 3 points)
     * @param centX   x coordinate of an optional point inside {@code
     *                graphXY} to use as the packing's alpha (centering)
     *                vertex; ignored unless {@code hasCent} is true
     * @param centY   y coordinate of that optional point; ignored unless
     *                {@code hasCent} is true
     * @param hasCent whether {@code centX}/{@code centY} should be used;
     *                pass false to let GOPack choose alpha automatically
     *                (matching {@code Packer::randomTri()}'s {@code
     *                cent=nullptr} default) -- also the fallback if the
     *                given point doesn't actually land inside {@code
     *                graphXY}
     * @return the generated complex and its raw (unpacked) layout. Throws
     *         {@link GOPackException} if generation fails (e.g. {@code
     *         graphXY} too short/degenerate, {@code intN}/{@code bdryN} too
     *         small to produce a usable complex). No riffle occurs, so there
     *         is no {@code maxPasses} parameter here -- it would have
     *         nothing to bound.
     */
    public static native RandomComplexResult computeRandomTri(
            int intN, int bdryN, double[] graphXY, double centX, double centY, boolean hasCent)
            throws GOPackException;

    /**
     * Generates a random Delaunay triangulation of the sphere and returns
     * its raw combinatorics, <b>not</b> a packing. This is the JNI bridge to
     * {@code gopack::Packer::randomTri(intN)} -- the sphere-only overload,
     * equivalent to {@code Packer::randomSphere(intN)} (see {@code
     * core/include/gopack/Packer.h}) -- with {@code setMode()}/{@code
     * riffle()} deliberately skipped, same as {@link #computeRandomTri}; see
     * the generator-group comment above that method.
     *
     * <p>No boundary/graph concept applies here (the sphere has no
     * boundary), so unlike {@link #computeRandomTri} there is no {@code
     * graphXY}/{@code cent} to pass. {@link RandomComplexResult#geometry} is
     * always Spherical ({@code +1}), and {@link RandomComplexResult#corners}
     * is always empty. <b>{@link RandomComplexResult#centersRe}/{@link
     * RandomComplexResult#centersIm} are meaningless here</b> (always {@code
     * 0}) -- the generator computes each point's (theta, phi) internally but
     * never carries it through to the returned complex (matching {@code
     * Packer::randomSphere}'s own behavior); only {@link
     * RandomComplexResult#flowers}/{@link RandomComplexResult#nodeCount}
     * matter for this method. {@link RandomComplexResult#radii} is likewise
     * just the usual uniform {@code 0.5} placeholder.
     *
     * @param intN number of points to generate (must be &gt;= 4)
     * @return the generated complex's raw combinatorics. Throws {@link
     *         GOPackException} if generation fails (e.g. {@code intN < 4}).
     */
    public static native RandomComplexResult computeRandomSphere(int intN)
            throws GOPackException;

    /**
     * Generates a random Delaunay triangulation of the rectangle
     * {@code [-aspect,aspect]x[-1,1]} and returns its raw combinatorics and
     * point layout, <b>not</b> a packing. This is the JNI bridge to
     * {@code gopack::Packer::randomRectangle(intN, aspect, bdryN)} (see
     * {@code core/include/gopack/Packer.h}); {@code riffle()} is
     * deliberately never called, same as {@link #computeRandomTri}; see the
     * generator-group comment above that method.
     *
     * <p>Unlike {@link #computeRandomTri}/{@link #computeRandomSphere}, the
     * native generator here does call {@code setMode(2, ...)} internally --
     * but only to identify and record the 4 actual corner vertices as
     * {@link RandomComplexResult#corners}; it still does not riffle, so the
     * result is not a packing. {@link RandomComplexResult#geometry} is
     * always Euclidean ({@code 0}). {@link RandomComplexResult#radii} is the
     * usual uniform {@code 0.5} placeholder; {@link
     * RandomComplexResult#centersRe}/{@link RandomComplexResult#centersIm}
     * hold the actual randomly-placed/Delaunay point positions.
     *
     * @param intN   number of interior points to generate (must be &gt;= 1)
     * @param aspect the rectangle's aspect ratio (half-width; height is
     *               fixed at {@code [-1,1]}); pass {@code 1.0} for a square
     * @param bdryN  number of boundary points to generate; pass {@code <= 0}
     *               to use the native default formula derived from {@code
     *               aspect} (matching {@code Packer::randomRectangle}'s own
     *               {@code bdryN=-1} default)
     * @return the generated complex and its raw (unpacked) layout. Throws
     *         {@link GOPackException} if generation fails (e.g. {@code
     *         intN < 1}).
     */
    public static native RandomComplexResult computeRandomRectangle(
            int intN, double aspect, int bdryN) throws GOPackException;

    /**
     * Generates a random Delaunay triangulation of the unit square and
     * returns its raw combinatorics and point layout, <b>not</b> a packing.
     * This is the JNI bridge to {@code gopack::Packer::randomSquare(n)} --
     * a thin wrapper that derives {@code intN}/{@code bdryN} from a single
     * total point count using its own interior/boundary split formula and
     * delegates to {@code Packer::randomRectangle(intN, 1.0, bdryN)} (see
     * {@code core/include/gopack/Packer.h}) -- exposed as its own bridge
     * rather than reimplementing that split formula in Java. Never packs;
     * see the generator-group comment above {@link #computeRandomTri}.
     *
     * <p>Same field semantics as {@link #computeRandomRectangle} with {@code
     * aspect=1.0}: {@link RandomComplexResult#geometry} is always Euclidean
     * ({@code 0}), {@link RandomComplexResult#corners} holds the 4 actual
     * corner vertices, {@link RandomComplexResult#radii} is the usual
     * uniform placeholder, and {@link RandomComplexResult#centersRe}/{@link
     * RandomComplexResult#centersIm} hold the actual point positions.
     *
     * @param n total number of points (interior + boundary) to generate
     * @return the generated complex and its raw (unpacked) layout. Throws
     *         {@link GOPackException} if generation fails (e.g. {@code n}
     *         too small to produce a usable complex).
     */
    public static native RandomComplexResult computeRandomSquare(int n)
            throws GOPackException;

    /**
     * Generates a random triangulation of the unit disc and computes its
     * maximal (hyperbolic) packing. This is the JNI bridge to {@code
     * gopack::Packer::randomDisc(N)} (see {@code core/include/gopack/Packer.h});
     * the CLI's {@code --random-disc} flag is the other way to reach the same
     * native code.
     *
     * <p>Like {@link #computeRandomTri}, the caller doesn't already know the
     * combinatorics -- the generator invents them from a Delaunay
     * triangulation of randomly placed points in the unit disc -- so this
     * returns a {@link RandomComplexResult} (flowers + radii + centers +
     * bookkeeping), not a bare {@code double[]} of radii. The result's {@link
     * RandomComplexResult#geometry} is always Hyperbolic ({@code -1}), and
     * {@link RandomComplexResult#corners} is always empty (the disc has no
     * polygonal corners the way a rectangle does).
     *
     * <p><b>Unlike {@link #computeRandomTri}/{@link #computeRandomSphere}/
     * {@link #computeRandomRectangle}/{@link #computeRandomSquare} above,
     * this method's result genuinely is a packing</b> ({@link
     * RandomComplexResult#radii}/{@link RandomComplexResult#centersRe}/
     * {@link RandomComplexResult#centersIm} come from a real {@code
     * riffle()}, not a placeholder) -- this backs the separate {@code
     * random_disc} command, which unlike {@code random_tri} deliberately
     * wants a max-packed disc, not a raw triangulation.
     *
     * @param n         total number of points (interior + boundary) to
     *                  generate; {@code Packer::randomDisc} splits this
     *                  automatically into roughly {@code sqrt(n)} boundary
     *                  points and the rest interior
     * @param maxPasses upper bound on riffle passes (pass &lt;= 0 to use the
     *                  native default of 200)
     * @return the generated complex and its computed hyperbolic packing.
     *         Throws {@link GOPackException} if generation fails (e.g. {@code
     *         n} too small to produce a usable complex) or riffle fails.
     */
    public static native RandomComplexResult computeRandomDisc(int n, int maxPasses)
            throws GOPackException;

    /** Returns the linked native library's version string, for diagnostics. */
    public static native String nativeVersion();
}
