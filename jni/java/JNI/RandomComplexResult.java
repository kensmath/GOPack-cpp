package JNI;

/**
 * The result of a native "geometrically random" triangulation generator
 * (e.g. {@link GOPackNative#computeRandomTri}), after it has been riffled to
 * a maximal (or polygonal) packing.
 *
 * <p>Unlike {@link GOPackNative#computeMaximalPacking} and
 * {@link GOPackNative#computeMaximalPackingFromComplex} -- where the caller
 * already supplies the combinatorics and only needs radii back -- a random
 * generator invents the combinatorics on the native side, so the caller
 * needs the full complex back: the flowers (to build its own vertex/complex
 * objects), and both the radii and centers of the computed packing (not just
 * radii), plus the bookkeeping fields (geometry, alpha, gamma, corners) that
 * describe how the packing was set up.
 *
 * <p>All per-vertex arrays ({@link #flowers}, {@link #radii},
 * {@link #centersRe}, {@link #centersIm}) are 1-indexed, length
 * {@code nodeCount+1} with index 0 unused/ignored -- the same convention
 * used throughout the native core and the rest of this JNI bridge, so a
 * caller never has to apply an off-by-one translation.
 */
public final class RandomComplexResult {

    /** Number of vertices in the generated complex. */
    public final int nodeCount;

    /**
     * v's petal list, 1-indexed ({@code flowers[0]} is an empty/unused
     * placeholder). Matches the *.p FLOWERS format convention: CLOSED (first
     * element == last element) iff v is an interior vertex, OPEN (first !=
     * last) iff v is a boundary vertex.
     */
    public final int[][] flowers;

    /** Euclidean radii after riffle(), 1-indexed, {@code radii[0]} unused. */
    public final double[] radii;

    /** Real part of each vertex's center after riffle(), 1-indexed, index 0 unused. */
    public final double[] centersRe;

    /** Imaginary part of each vertex's center after riffle(), 1-indexed, index 0 unused. */
    public final double[] centersIm;

    /**
     * 0 = Euclidean, -1 = Hyperbolic, +1 = Spherical -- matches
     * {@code gopack::Geometry}'s underlying values and the *.p file's
     * GEOMETRY: field, same convention already used by
     * {@link GOPackNative#computeMaximalPackingFromComplex}'s geometry
     * parameter.
     */
    public final int geometry;

    /** The packing's alpha (centering) vertex, chosen by the generator. */
    public final int alpha;

    /** The packing's gamma vertex (if any; 0 if unset), chosen by the generator. */
    public final int gamma;

    /**
     * Corner vertex numbers, in counterclockwise order, for a polygonal
     * (rectangle/square) generator; empty for generators that stay in
     * max-pack mode (e.g. {@link GOPackNative#computeRandomTri}).
     */
    public final int[] corners;

    public RandomComplexResult(int nodeCount, int[][] flowers, double[] radii,
            double[] centersRe, double[] centersIm, int geometry, int alpha, int gamma,
            int[] corners) {
        this.nodeCount = nodeCount;
        this.flowers = flowers;
        this.radii = radii;
        this.centersRe = centersRe;
        this.centersIm = centersIm;
        this.geometry = geometry;
        this.alpha = alpha;
        this.gamma = gamma;
        this.corners = corners;
    }
}
