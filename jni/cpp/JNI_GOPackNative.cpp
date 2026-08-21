// JNI implementation backing JNI.GOPackNative (Java package "JNI", renamed
// from org.kensmath.gopack 8/2026 -- see GOPackNative.java's doc comment).
//
// This file is a thin marshalling layer only: it converts JNI arguments to
// C++ types, calls into gopack_core's Packer class, and converts results/
// exceptions back. All actual packing logic lives in core/.
#include <jni.h>

#include <string>

#include "gopack/Packer.h"

namespace {

void throwGOPackException(JNIEnv* env, const std::string& message) {
    jclass exClass = env->FindClass("JNI/GOPackException");
    if (exClass != nullptr) {
        env->ThrowNew(exClass, message.c_str());
    }
}

std::string jstringToString(JNIEnv* env, jstring js) {
    const char* chars = env->GetStringUTFChars(js, nullptr);
    std::string s(chars);
    env->ReleaseStringUTFChars(js, chars);
    return s;
}

#ifdef GOPACK_HAVE_RANDOM_GEN
// Marshals a fully-riffled Packer (after a random generator + riffle() call)
// into a new JNI.RandomComplexResult, carrying everything a
// caller needs to reconstruct the complex on the Java side (flowers) as well
// as the resulting packing (radii/centers) -- unlike
// computeMaximalPacking[FromComplex] above, a random generator's caller
// doesn't already know the combinatorics, so a bare double[] of radii isn't
// enough here (see RandomComplexResult.java for the full field rationale).
// Returns nullptr, with a pending Java exception already thrown via
// throwGOPackException, on failure (matching this file's other helpers).
jobject buildRandomComplexResult(JNIEnv* env, const gopack::Packer& packer) {
    jclass resultClass = env->FindClass("JNI/RandomComplexResult");
    if (resultClass == nullptr) {
        throwGOPackException(env,
            "GOPack native: JNI.RandomComplexResult class not found -- is "
            "it on the classpath?");
        return nullptr;
    }
    jmethodID ctor = env->GetMethodID(resultClass, "<init>", "(I[[I[D[D[DIII[I)V");
    if (ctor == nullptr) {
        throwGOPackException(env,
            "GOPack native: RandomComplexResult constructor not found (signature mismatch?)");
        return nullptr;
    }

    const jsize n = static_cast<jsize>(packer.nodeCount) + 1;

    // flowers: 1-indexed jagged int[][], flowers[0] left as an empty (not
    // null) int[] placeholder so a naive 0..nodeCount Java-side loop can't
    // NPE on it.
    jclass intArrayClass = env->FindClass("[I");
    jobjectArray flowersArr = env->NewObjectArray(n, intArrayClass, nullptr);
    for (jsize v = 0; v < n; ++v) {
        static const std::vector<gopack::Index> kEmpty;
        const std::vector<gopack::Index>& petals =
            (v == 0) ? kEmpty : packer.flowers[static_cast<size_t>(v)];
        jintArray row = env->NewIntArray(static_cast<jsize>(petals.size()));
        if (!petals.empty()) {
            env->SetIntArrayRegion(row, 0, static_cast<jsize>(petals.size()), petals.data());
        }
        env->SetObjectArrayElement(flowersArr, v, row);
        env->DeleteLocalRef(row);
    }

    jdoubleArray radiiArr = env->NewDoubleArray(n);
    env->SetDoubleArrayRegion(radiiArr, 0, n, packer.radii.data());

    // packer.centers is std::vector<std::complex<double>> -- split into
    // parallel re/im double[] arrays, since JNI has no complex type.
    std::vector<jdouble> re(static_cast<size_t>(n)), im(static_cast<size_t>(n));
    for (jsize v = 0; v < n; ++v) {
        re[static_cast<size_t>(v)] = packer.centers[static_cast<size_t>(v)].real();
        im[static_cast<size_t>(v)] = packer.centers[static_cast<size_t>(v)].imag();
    }
    jdoubleArray centersReArr = env->NewDoubleArray(n);
    env->SetDoubleArrayRegion(centersReArr, 0, n, re.data());
    jdoubleArray centersImArr = env->NewDoubleArray(n);
    env->SetDoubleArrayRegion(centersImArr, 0, n, im.data());

    jintArray cornersArr = env->NewIntArray(static_cast<jsize>(packer.corners.size()));
    if (!packer.corners.empty()) {
        env->SetIntArrayRegion(cornersArr, 0, static_cast<jsize>(packer.corners.size()),
                                packer.corners.data());
    }

    return env->NewObject(resultClass, ctor, static_cast<jint>(packer.nodeCount), flowersArr,
                           radiiArr, centersReArr, centersImArr,
                           static_cast<jint>(packer.hes), static_cast<jint>(packer.alpha),
                           static_cast<jint>(packer.gamma), cornersArr);
}
#endif // GOPACK_HAVE_RANDOM_GEN

} // namespace

extern "C" {

// double[] computeMaximalPacking(String inputPath, int geometryHint, double
// tolerance (currently informational only -- see NOTE below), int
// maxPasses)
//
// Returns the euclidean radii for every vertex, as an array of length
// nodeCount+1 with radii[v] holding vertex v's radius for v = 1..nodeCount
// (matching GOPack/CirclePack's own 1-indexed vertex numbering, and the
// same convention core/include/gopack/Packer.h uses throughout); radii[0]
// is unused/unspecified. This is a deliberate choice over a "natural" Java
// 0-indexed array with radii[i] = vertex (i+1)'s radius: keeping vertex v
// at index v everywhere -- MATLAB, the C++ core, and now this bridge --
// means a caller never has to remember or apply an off-by-one translation.
// NOTE: 'tolerance' is accepted for forward compatibility with the
// eventual configurable stopping criterion, but the current port uses
// GOPack's fixed 0.01 visual-error cutoff (continueRiffle.m's 'cutval'); a
// future pass can thread a real tolerance through once that's wired up.
JNIEXPORT jdoubleArray JNICALL
Java_JNI_GOPackNative_computeMaximalPacking(
    JNIEnv* env, jclass /*clazz*/, jstring inputPath, jint /*geometryHint*/,
    jdouble /*tolerance*/, jint maxPasses) {

    try {
        const std::string path = jstringToString(env, inputPath);

        gopack::Packer packer;
        if (packer.readpack(path) <= 0) {
            throwGOPackException(env, "GOPack native: failed to read input file: " + path);
            return nullptr;
        }
        if (packer.setMode(1) < 0) {
            throwGOPackException(env, "GOPack native: failed to set max-pack mode");
            return nullptr;
        }

        // 200, not the old 20, matching cli/main.cpp's own default -- 20 was
        // found to under-converge small/irregular random triangulations (see
        // that file's comment); continueRiffle() exits early on convergence
        // regardless, so a higher cap costs nothing for inputs that finish
        // sooner.
        gopack::RiffleResult result = packer.riffle(maxPasses > 0 ? maxPasses : 200);
        if (result.cycles < 0) {
            throwGOPackException(env, "GOPack native: riffle failed");
            return nullptr;
        }

        // packer.radii is 1-indexed (size nodeCount+1, index 0 unused);
        // return it as-is (same size, same indexing) rather than shifting
        // to a 0-indexed Java array, so radii[v] is always vertex v's
        // radius on both sides of the JNI boundary.
        const jsize n = static_cast<jsize>(packer.nodeCount) + 1;
        jdoubleArray out = env->NewDoubleArray(n);
        if (out == nullptr) {
            throwGOPackException(env, "Failed to allocate result array");
            return nullptr;
        }
        env->SetDoubleArrayRegion(out, 0, n, packer.radii.data());
        return out;
    } catch (const std::exception& e) {
        throwGOPackException(env, e.what());
        return nullptr;
    }
}

// double[][] computeMaximalPackingFromComplex(int nodeCount, int[][] flowers,
// int geometry, double tolerance, int maxPasses)
//
// The in-memory counterpart to computeMaximalPacking above: takes a
// triangulation already held in memory by the caller (e.g. CirclePack's own
// per-vertex flower data) instead of a file path, so the combinatorics never
// have to round-trip through the *.p text format just to be handed straight
// back to this native call. This is the entry point Packer::loadComplex()
// (core/src/PackerIO.cpp) exists for; see that method's doc comment for the
// full semantics.
//
// flowers: length nodeCount+1, 1-indexed (flowers[0] is ignored/unused).
// flowers[v] is v's petal list in the *.p FLOWERS format convention: CLOSED
// (first element == last element) iff v is an interior vertex, OPEN
// (first != last) iff v is a boundary vertex.
// geometry: 0 = Euclidean, -1 = Hyperbolic, +1 = Spherical (matches
// gopack::Geometry's own underlying values, and GOPack's GEOMETRY: file
// field).
//
// Returns a 3-row double[][], each row of length nodeCount+1 (1-indexed,
// index 0 unused), matching GOPackNative.java's documented contract:
//   result[0] -- radii
//   result[1] -- center real parts
//   result[2] -- center imaginary parts
// IMPORTANT: this MUST stay a jobjectArray of 3 jdoubleArrays, matching the
// Java side's "public static native double[][] ..." declaration exactly.
// JNI resolves/links native methods purely by name+argument-descriptor --
// it does NOT check the native function's actual return type against the
// declared Java return type. A JNIEXPORT here that instead returned a bare
// jdoubleArray (as an earlier version of this function mistakenly did)
// would still link and run without any JNI-level error, but the JVM
// interpreter -- trusting the *Java-declared* double[][] type -- would
// execute result[0]/result[1]/result[2] as aaload (reference-array) reads
// against what is actually raw double payload data, reinterpreting radius/
// center bits as compressed-oop pointers. That is a real regression this
// project hit once already (see git history around the HypPacker.maxPack
// EXCEPTION_ACCESS_VIOLATION crash, N=10000): the corrupted "pointer" was
// literally the high 32 bits of an ordinary ~0.5 radius value. Do not
// "simplify" this back to a flat jdoubleArray.
JNIEXPORT jobjectArray JNICALL
Java_JNI_GOPackNative_computeMaximalPackingFromComplex(
    JNIEnv* env, jclass /*clazz*/, jint nodeCount, jobjectArray flowers, jint geometry,
    jdouble /*tolerance*/, jint maxPasses) {

    try {
        if (nodeCount <= 0) {
            throwGOPackException(env, "GOPack native: nodeCount must be positive");
            return nullptr;
        }
        const jsize expectedLen = static_cast<jsize>(nodeCount) + 1;
        if (env->GetArrayLength(flowers) != expectedLen) {
            throwGOPackException(env,
                "GOPack native: flowers.length must be nodeCount+1 (index 0 unused)");
            return nullptr;
        }

        std::vector<std::vector<gopack::Index>> flowersIn(static_cast<size_t>(expectedLen));
        for (jsize v = 1; v < expectedLen; ++v) {
            jobject rowObj = env->GetObjectArrayElement(flowers, v);
            if (rowObj == nullptr) {
                throwGOPackException(env,
                    "GOPack native: flowers[v] must not be null for v=1..nodeCount");
                return nullptr;
            }
            jintArray row = static_cast<jintArray>(rowObj);
            jsize rowLen = env->GetArrayLength(row);
            jint* rowData = env->GetIntArrayElements(row, nullptr);
            std::vector<gopack::Index> flower(rowData, rowData + rowLen);
            env->ReleaseIntArrayElements(row, rowData, JNI_ABORT);
            env->DeleteLocalRef(rowObj);
            flowersIn[static_cast<size_t>(v)] = std::move(flower);
        }

        gopack::Packer packer;
        gopack::Index loaded = packer.loadComplex(static_cast<gopack::Index>(nodeCount),
                                                    flowersIn,
                                                    static_cast<gopack::Geometry>(geometry));
        if (loaded <= 0) {
            throwGOPackException(env,
                "GOPack native: loadComplex failed (invalid/malformed complex -- see stderr)");
            return nullptr;
        }
        if (packer.setMode(1) < 0) {
            throwGOPackException(env, "GOPack native: failed to set max-pack mode");
            return nullptr;
        }

        // 200, not the old 20, matching cli/main.cpp's own default -- 20 was
        // found to under-converge small/irregular random triangulations (see
        // that file's comment); continueRiffle() exits early on convergence
        // regardless, so a higher cap costs nothing for inputs that finish
        // sooner.
        gopack::RiffleResult result = packer.riffle(maxPasses > 0 ? maxPasses : 200);
        if (result.cycles < 0) {
            throwGOPackException(env, "GOPack native: riffle failed");
            return nullptr;
        }

        // packer.centers is std::vector<std::complex<double>> -- JNI has no
        // complex type, so split into parallel re/im double[] rows (same
        // approach buildRandomComplexResult() uses below for the random
        // generators' RandomComplexResult.centersRe/centersIm fields).
        const jsize n = static_cast<jsize>(packer.nodeCount) + 1;

        jdoubleArray radiiArr = env->NewDoubleArray(n);
        if (radiiArr == nullptr) {
            throwGOPackException(env, "Failed to allocate radii result array");
            return nullptr;
        }
        env->SetDoubleArrayRegion(radiiArr, 0, n, packer.radii.data());

        std::vector<jdouble> re(static_cast<size_t>(n)), im(static_cast<size_t>(n));
        for (jsize v = 0; v < n; ++v) {
            re[static_cast<size_t>(v)] = packer.centers[static_cast<size_t>(v)].real();
            im[static_cast<size_t>(v)] = packer.centers[static_cast<size_t>(v)].imag();
        }
        jdoubleArray centerXArr = env->NewDoubleArray(n);
        jdoubleArray centerYArr = env->NewDoubleArray(n);
        if (centerXArr == nullptr || centerYArr == nullptr) {
            throwGOPackException(env, "Failed to allocate centers result array");
            return nullptr;
        }
        env->SetDoubleArrayRegion(centerXArr, 0, n, re.data());
        env->SetDoubleArrayRegion(centerYArr, 0, n, im.data());

        jclass doubleArrayClass = env->FindClass("[D");
        if (doubleArrayClass == nullptr) {
            throwGOPackException(env, "GOPack native: [D class not found");
            return nullptr;
        }
        jobjectArray out = env->NewObjectArray(3, doubleArrayClass, nullptr);
        if (out == nullptr) {
            throwGOPackException(env, "Failed to allocate result array");
            return nullptr;
        }
        env->SetObjectArrayElement(out, 0, radiiArr);
        env->SetObjectArrayElement(out, 1, centerXArr);
        env->SetObjectArrayElement(out, 2, centerYArr);
        return out;
    } catch (const std::exception& e) {
        throwGOPackException(env, e.what());
        return nullptr;
    }
}

// double[][] computePolygonalPackingFromComplex(int nodeCount, int[][] flowers,
// int geometry, int[] corners, double[] angles, int maxPasses)
//
// Polygonal/rectangle packing (mode 2) counterpart to
// computeMaximalPackingFromComplex above: loads a triangulation already held
// in memory (same flowers/geometry convention as that method -- see its own
// comment), then instead of maximal-packing mode, sets polygonal mode via
// Packer::setMode(2, corners, angles) before riffling. See setMode()'s doc
// comment in Packer.h for the authoritative corner/angle semantics; the
// short version:
//
// corners: 1-indexed boundary vertices to use as polygon corners, in
// counterclockwise order. Pass null or a 0-length array to let GOPack infer/
// choose corners automatically (mirrors Packer::setMode's own crns={}
// default -- "figure out the corners for me").
//
// angles: target corner angles as MULTIPLES OF PI (the same convention
// CirclePack's own set_aim command uses -- e.g. pass 0.5 for a right angle,
// NOT Math.PI/2.0 -- chosen for consistency with that command and because
// it's the easier convention for interactive/user input), matching corners
// in count and order. This bridge converts to radians internally before
// calling Packer::setMode(), which itself works in radians. Pass null or a
// 0-length array for automatic equal angles (e.g. 0.5 each for a 4-corner
// rectangle -- setMode's own default). Must be empty whenever corners is
// empty; passing angles without corners is rejected, since there'd be no
// way to say which angle belongs to which (as-yet-unchosen) corner. The
// angles must be consistent with the polygon's turning-angle requirement
// (sum of (1 - angle) over all corners == 2, i.e. in radians sum of
// (pi - angle) == 2*pi); setMode() itself validates this (in radians, after
// this bridge's conversion) and this bridge surfaces that failure as a
// GOPackException rather than silently producing a bad packing.
//
// NOTE: entering mode 2 always resets the packing's internal geometry to
// Euclidean, regardless of the 'geometry' passed in here (a deliberate
// GOPack behavior -- see setMode()'s own .cpp comment) -- polygonal
// packings are inherently planar.
//
// Returns a 4-row double[][], each row of length nodeCount+1 (1-indexed,
// index 0 unused):
//   result[0] -- radii
//   result[1] -- center real parts
//   result[2] -- center imaginary parts
//   result[3] -- the actual corner vertices used, as doubles (exact for the
//                small integers involved) -- this is what you gave in
//                'corners' when non-empty, or GOPack's own automatic choice
//                when 'corners' was null/empty; length == number of
//                corners actually used (NOT nodeCount+1 like the other three
//                rows -- check its length rather than assuming 4).
JNIEXPORT jobjectArray JNICALL
Java_JNI_GOPackNative_computePolygonalPackingFromComplex(
    JNIEnv* env, jclass /*clazz*/, jint nodeCount, jobjectArray flowers, jint geometry,
    jintArray corners, jdoubleArray angles, jint maxPasses) {

    try {
        if (nodeCount <= 0) {
            throwGOPackException(env, "GOPack native: nodeCount must be positive");
            return nullptr;
        }
        const jsize expectedLen = static_cast<jsize>(nodeCount) + 1;
        if (env->GetArrayLength(flowers) != expectedLen) {
            throwGOPackException(env,
                "GOPack native: flowers.length must be nodeCount+1 (index 0 unused)");
            return nullptr;
        }

        std::vector<std::vector<gopack::Index>> flowersIn(static_cast<size_t>(expectedLen));
        for (jsize v = 1; v < expectedLen; ++v) {
            jobject rowObj = env->GetObjectArrayElement(flowers, v);
            if (rowObj == nullptr) {
                throwGOPackException(env,
                    "GOPack native: flowers[v] must not be null for v=1..nodeCount");
                return nullptr;
            }
            jintArray row = static_cast<jintArray>(rowObj);
            jsize rowLen = env->GetArrayLength(row);
            jint* rowData = env->GetIntArrayElements(row, nullptr);
            std::vector<gopack::Index> flower(rowData, rowData + rowLen);
            env->ReleaseIntArrayElements(row, rowData, JNI_ABORT);
            env->DeleteLocalRef(rowObj);
            flowersIn[static_cast<size_t>(v)] = std::move(flower);
        }

        gopack::Packer packer;
        gopack::Index loaded = packer.loadComplex(static_cast<gopack::Index>(nodeCount),
                                                    flowersIn,
                                                    static_cast<gopack::Geometry>(geometry));
        if (loaded <= 0) {
            throwGOPackException(env,
                "GOPack native: loadComplex failed (invalid/malformed complex -- see stderr)");
            return nullptr;
        }

        // corners: 1-indexed boundary vertices, cclw order; null/empty ->
        // let GOPack infer them (Packer::setMode's own crns={} default).
        std::vector<gopack::Index> cornersIn;
        if (corners != nullptr && env->GetArrayLength(corners) > 0) {
            const jsize cornersLen = env->GetArrayLength(corners);
            jint* cornersData = env->GetIntArrayElements(corners, nullptr);
            cornersIn.assign(cornersData, cornersData + cornersLen);
            env->ReleaseIntArrayElements(corners, cornersData, JNI_ABORT);
        }

        // angles: target corner angles as multiples of pi (CirclePack's own
        // set_aim convention -- see this function's doc comment above),
        // matching corners in count and order; null/empty -> automatic
        // equal angles. Convert to radians here, since Packer::setMode()
        // (like the MATLAB setMode.m it's ported from) works in radians.
        std::vector<gopack::Scalar> anglesIn;
        if (angles != nullptr && env->GetArrayLength(angles) > 0) {
            const jsize anglesLen = env->GetArrayLength(angles);
            jdouble* anglesData = env->GetDoubleArrayElements(angles, nullptr);
            anglesIn.resize(static_cast<size_t>(anglesLen));
            for (jsize i = 0; i < anglesLen; ++i) {
                anglesIn[static_cast<size_t>(i)] = anglesData[i] * gopack::kPi;
            }
            env->ReleaseDoubleArrayElements(angles, anglesData, JNI_ABORT);
        }
        if (!anglesIn.empty() && cornersIn.empty()) {
            throwGOPackException(env,
                "GOPack native: angles given without corners -- corners must be given "
                "explicitly to pair with angles");
            return nullptr;
        }
        if (!anglesIn.empty() && anglesIn.size() != cornersIn.size()) {
            throwGOPackException(env,
                "GOPack native: angles.length must match corners.length (or both be empty)");
            return nullptr;
        }

        if (packer.setMode(2, cornersIn, anglesIn) < 0) {
            throwGOPackException(env,
                "GOPack native: setMode(2, corners, angles) failed -- see stderr for the "
                "specific diagnostic (e.g. a given corner isn't a boundary vertex, or corner "
                "angles aren't consistent with the polygon's turning-angle requirement)");
            return nullptr;
        }

        gopack::RiffleResult result = packer.riffle(maxPasses > 0 ? maxPasses : 200);
        if (result.cycles < 0) {
            throwGOPackException(env, "GOPack native: riffle failed");
            return nullptr;
        }

        // Same radii/centerX/centerY marshalling as
        // computeMaximalPackingFromComplex above, plus a 4th row for the
        // corner vertices actually used (see this function's doc comment).
        const jsize n = static_cast<jsize>(packer.nodeCount) + 1;

        jdoubleArray radiiArr = env->NewDoubleArray(n);
        if (radiiArr == nullptr) {
            throwGOPackException(env, "Failed to allocate radii result array");
            return nullptr;
        }
        env->SetDoubleArrayRegion(radiiArr, 0, n, packer.radii.data());

        std::vector<jdouble> re(static_cast<size_t>(n)), im(static_cast<size_t>(n));
        for (jsize v = 0; v < n; ++v) {
            re[static_cast<size_t>(v)] = packer.centers[static_cast<size_t>(v)].real();
            im[static_cast<size_t>(v)] = packer.centers[static_cast<size_t>(v)].imag();
        }
        jdoubleArray centerXArr = env->NewDoubleArray(n);
        jdoubleArray centerYArr = env->NewDoubleArray(n);
        if (centerXArr == nullptr || centerYArr == nullptr) {
            throwGOPackException(env, "Failed to allocate centers result array");
            return nullptr;
        }
        env->SetDoubleArrayRegion(centerXArr, 0, n, re.data());
        env->SetDoubleArrayRegion(centerYArr, 0, n, im.data());

        const jsize cornerCount = static_cast<jsize>(packer.corners.size());
        std::vector<jdouble> cornersOut(static_cast<size_t>(cornerCount));
        for (jsize i = 0; i < cornerCount; ++i) {
            cornersOut[static_cast<size_t>(i)] = static_cast<jdouble>(packer.corners[static_cast<size_t>(i)]);
        }
        jdoubleArray cornersArr = env->NewDoubleArray(cornerCount);
        if (cornersArr == nullptr) {
            throwGOPackException(env, "Failed to allocate corners result array");
            return nullptr;
        }
        if (cornerCount > 0) {
            env->SetDoubleArrayRegion(cornersArr, 0, cornerCount, cornersOut.data());
        }

        jclass doubleArrayClass = env->FindClass("[D");
        if (doubleArrayClass == nullptr) {
            throwGOPackException(env, "GOPack native: [D class not found");
            return nullptr;
        }
        jobjectArray out = env->NewObjectArray(4, doubleArrayClass, nullptr);
        if (out == nullptr) {
            throwGOPackException(env, "Failed to allocate result array");
            return nullptr;
        }
        env->SetObjectArrayElement(out, 0, radiiArr);
        env->SetObjectArrayElement(out, 1, centerXArr);
        env->SetObjectArrayElement(out, 2, centerYArr);
        env->SetObjectArrayElement(out, 3, cornersArr);
        return out;
    } catch (const std::exception& e) {
        throwGOPackException(env, e.what());
        return nullptr;
    }
}

#ifdef GOPACK_HAVE_RANDOM_GEN
// RandomComplexResult computeRandomTri(int intN, int bdryN, double[] graphXY,
//     double centX, double centY, boolean hasCent, int maxPasses)
//
// Bridges Packer::randomTri(intN, bdryN, graph, cent) -- a random
// triangulation of an ARBITRARY closed polygonal region, not just a fixed
// disc/square/rectangle -- so CirclePack can hand it a user-drawn boundary
// and get back a full ready-to-render complex. See
// RandomComplexResult.java and this file's computeRandomTri Javadoc
// counterpart in GOPackNative.java for the full field semantics.
//
// graphXY: flat x,y coordinate list for the closed boundary polygon (do not
// repeat the first point at the end); must have an even length of at least
// 6 (i.e. at least 3 points).
JNIEXPORT jobject JNICALL
Java_JNI_GOPackNative_computeRandomTri(
    JNIEnv* env, jclass /*clazz*/, jint intN, jint bdryN, jdoubleArray graphXY, jdouble centX,
    jdouble centY, jboolean hasCent, jint maxPasses) {

    try {
        const jsize flatLen = env->GetArrayLength(graphXY);
        if (flatLen < 6 || (flatLen % 2) != 0) {
            throwGOPackException(env,
                "GOPack native: graphXY must hold at least 3 (x,y) points and have an even "
                "length");
            return nullptr;
        }
        jdouble* flat = env->GetDoubleArrayElements(graphXY, nullptr);
        std::vector<gopack::Complex> graph;
        graph.reserve(static_cast<size_t>(flatLen) / 2);
        for (jsize i = 0; i + 1 < flatLen; i += 2) {
            graph.emplace_back(flat[i], flat[i + 1]);
        }
        env->ReleaseDoubleArrayElements(graphXY, flat, JNI_ABORT);

        gopack::Complex cent(centX, centY);
        gopack::Packer packer = gopack::Packer::randomTri(
            static_cast<gopack::Index>(intN), static_cast<gopack::Index>(bdryN), graph,
            hasCent ? &cent : nullptr);
        if (packer.nodeCount <= 0) {
            throwGOPackException(env,
                "GOPack native: randomTri failed to produce a usable complex (bad graph, or "
                "intN/bdryN too small -- see stderr)");
            return nullptr;
        }
        // randomTri() leaves mode at its Packer()-default of 1 (max-pack)
        // but does NOT itself call setMode() -- unlike randomRectangle()/
        // randomSquare(), which already call setMode(2, ...) internally --
        // so vAims (riffle()'s target angle sums) still needs to be
        // populated here, exactly as cli/main.cpp's own --random-tri
        // handling and the two bridges above already do for their sources.
        if (packer.setMode(1) < 0) {
            throwGOPackException(env, "GOPack native: failed to set max-pack mode");
            return nullptr;
        }

        gopack::RiffleResult result = packer.riffle(maxPasses > 0 ? maxPasses : 200);
        if (result.cycles < 0) {
            throwGOPackException(env, "GOPack native: riffle failed");
            return nullptr;
        }

        return buildRandomComplexResult(env, packer);
    } catch (const std::exception& e) {
        throwGOPackException(env, e.what());
        return nullptr;
    }
}
#endif // GOPACK_HAVE_RANDOM_GEN

#ifdef GOPACK_HAVE_RANDOM_GEN
// RandomComplexResult computeRandomTriLayout(int intN, int bdryN, double[]
//     graphXY, double centX, double centY, boolean hasCent)
//
// Raw-layout counterpart to computeRandomTri() above: same generator
// (Packer::randomTri(intN, bdryN, graph, cent) -- a Delaunay triangulation
// of randomly placed points inside an arbitrary closed polygonal region),
// but deliberately skips setMode(1)/riffle() so the result is the
// triangulation's own raw post-Delaunay layout (combinatorics + the actual
// Euclidean positions the random points were placed at) instead of that
// triangulation's intrinsic maximal packing.
//
// Why this exists as a separate method rather than a flag on
// computeRandomTri: for a disc-topology complex, the maximal packing is,
// by construction/uniformization, the canonical packing that fills the
// unit disc -- completely independent of graphXY's actual Euclidean shape.
// A caller that wants the result to visually resemble the input region
// (e.g. CirclePack's own pure-Java RandomTriangulation/Triangulation path,
// which this bridge is meant to eventually replace for large complexes)
// needs the raw layout, not a repacked one; the two are not reachable from
// each other by any post-hoc transform of the riffled result, since
// riffle() already discards the original positions. See
// HANDOFFrandomtrinorepack.md for the full request this implements.
//
// graphXY: same convention as computeRandomTri -- flat x,y coordinate list
// for the closed boundary polygon (do not repeat the first point at the
// end); must have an even length of at least 6 (i.e. at least 3 points).
//
// No maxPasses parameter: unlike computeRandomTri, nothing here ever
// riffles, so a pass-count bound would be silently ignored -- omitted
// entirely rather than kept as a confusing no-op argument.
JNIEXPORT jobject JNICALL
Java_JNI_GOPackNative_computeRandomTriLayout(
    JNIEnv* env, jclass /*clazz*/, jint intN, jint bdryN, jdoubleArray graphXY, jdouble centX,
    jdouble centY, jboolean hasCent) {

    try {
        const jsize flatLen = env->GetArrayLength(graphXY);
        if (flatLen < 6 || (flatLen % 2) != 0) {
            throwGOPackException(env,
                "GOPack native: graphXY must hold at least 3 (x,y) points and have an even "
                "length");
            return nullptr;
        }
        jdouble* flat = env->GetDoubleArrayElements(graphXY, nullptr);
        std::vector<gopack::Complex> graph;
        graph.reserve(static_cast<size_t>(flatLen) / 2);
        for (jsize i = 0; i + 1 < flatLen; i += 2) {
            graph.emplace_back(flat[i], flat[i + 1]);
        }
        env->ReleaseDoubleArrayElements(graphXY, flat, JNI_ABORT);

        gopack::Complex cent(centX, centY);
        gopack::Packer packer = gopack::Packer::randomTri(
            static_cast<gopack::Index>(intN), static_cast<gopack::Index>(bdryN), graph,
            hasCent ? &cent : nullptr);
        if (packer.nodeCount <= 0) {
            throwGOPackException(env,
                "GOPack native: randomTri failed to produce a usable complex (bad graph, or "
                "intN/bdryN too small -- see stderr)");
            return nullptr;
        }

        // Deliberately no setMode()/riffle() here -- see this function's
        // header comment. randomTri() itself already leaves the packer with
        // hes = Euclidean and radii/centers populated by parseTriangles()
        // (radii: a uniform 0.5 placeholder with no packing meaning;
        // centers: the actual randomly-placed/Delaunay point positions) --
        // buildRandomComplexResult() reads those fields directly and has no
        // riffle()-completion assumption of its own, so it works unchanged
        // on this non-riffled packer.
        return buildRandomComplexResult(env, packer);
    } catch (const std::exception& e) {
        throwGOPackException(env, e.what());
        return nullptr;
    }
}
#endif // GOPACK_HAVE_RANDOM_GEN

#ifdef GOPACK_HAVE_RANDOM_GEN
// RandomComplexResult computeRandomDisc(int n, int maxPasses)
//
// Bridges Packer::randomDisc(N) -- a random triangulation of the unit disc,
// set up for hyperbolic maximal packing -- so CirclePack can generate a
// fresh random disc packing without round-tripping through a *.p file. See
// RandomComplexResult.java and this file's computeRandomDisc Javadoc
// counterpart in GOPackNative.java for the full field semantics.
JNIEXPORT jobject JNICALL
Java_JNI_GOPackNative_computeRandomDisc(
    JNIEnv* env, jclass /*clazz*/, jint n, jint maxPasses) {

    try {
        gopack::Packer packer = gopack::Packer::randomDisc(static_cast<gopack::Index>(n));
        if (packer.nodeCount <= 0) {
            throwGOPackException(env,
                "GOPack native: randomDisc failed to produce a usable complex (n too small? "
                "-- see stderr)");
            return nullptr;
        }
        // randomDisc() sets hes=Hyperbolic and mode=1 directly on the struct
        // but, like randomTri(), does NOT itself call setMode() -- so vAims
        // (riffle()'s target angle sums) still needs to be populated here.
        // (randomSquare()/randomRectangle() are the ones that already call
        // setMode(2, ...) internally -- calling setMode(1) on those would be
        // wrong; randomDisc() is not one of those, confirmed directly in
        // core/src/PackerRandom.cpp.)
        if (packer.setMode(1) < 0) {
            throwGOPackException(env, "GOPack native: failed to set max-pack mode");
            return nullptr;
        }

        gopack::RiffleResult result = packer.riffle(maxPasses > 0 ? maxPasses : 200);
        if (result.cycles < 0) {
            throwGOPackException(env, "GOPack native: riffle failed");
            return nullptr;
        }

        return buildRandomComplexResult(env, packer);
    } catch (const std::exception& e) {
        throwGOPackException(env, e.what());
        return nullptr;
    }
}
#endif // GOPACK_HAVE_RANDOM_GEN

JNIEXPORT jstring JNICALL
Java_JNI_GOPackNative_nativeVersion(JNIEnv* env, jclass /*clazz*/) {
    return env->NewStringUTF(
        "gopack-cpp 0.1.0 (max-pack mode ported; polygonal/rectangle modes pending)");
}

} // extern "C"
