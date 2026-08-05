package org.kensmath.gopack;

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
     * @param maxPasses    upper bound on riffle passes (GOPack default is 20)
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

    /** Returns the linked native library's version string, for diagnostics. */
    public static native String nativeVersion();
}
