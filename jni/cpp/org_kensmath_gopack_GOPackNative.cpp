// JNI implementation backing org.kensmath.gopack.GOPackNative.
//
// This file is a thin marshalling layer only: it converts JNI arguments to
// C++ types, calls into gopack_core's Packer class, and converts results/
// exceptions back. All actual packing logic lives in core/.
#include <jni.h>

#include <string>

#include "gopack/Packer.h"

namespace {

void throwGOPackException(JNIEnv* env, const std::string& message) {
    jclass exClass = env->FindClass("org/kensmath/gopack/GOPackException");
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
// into a new org.kensmath.gopack.RandomComplexResult, carrying everything a
// caller needs to reconstruct the complex on the Java side (flowers) as well
// as the resulting packing (radii/centers) -- unlike
// computeMaximalPacking[FromComplex] above, a random generator's caller
// doesn't already know the combinatorics, so a bare double[] of radii isn't
// enough here (see RandomComplexResult.java for the full field rationale).
// Returns nullptr, with a pending Java exception already thrown via
// throwGOPackException, on failure (matching this file's other helpers).
jobject buildRandomComplexResult(JNIEnv* env, const gopack::Packer& packer) {
    jclass resultClass = env->FindClass("org/kensmath/gopack/RandomComplexResult");
    if (resultClass == nullptr) {
        throwGOPackException(env,
            "GOPack native: org.kensmath.gopack.RandomComplexResult class not found -- is "
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
Java_org_kensmath_gopack_GOPackNative_computeMaximalPacking(
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

// double[] computeMaximalPackingFromComplex(int nodeCount, int[][] flowers,
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
// Returns radii in the same nodeCount+1-length, 1-indexed convention as
// computeMaximalPacking.
JNIEXPORT jdoubleArray JNICALL
Java_org_kensmath_gopack_GOPackNative_computeMaximalPackingFromComplex(
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
Java_org_kensmath_gopack_GOPackNative_computeRandomTri(
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
// RandomComplexResult computeRandomDisc(int n, int maxPasses)
//
// Bridges Packer::randomDisc(N) -- a random triangulation of the unit disc,
// set up for hyperbolic maximal packing -- so CirclePack can generate a
// fresh random disc packing without round-tripping through a *.p file. See
// RandomComplexResult.java and this file's computeRandomDisc Javadoc
// counterpart in GOPackNative.java for the full field semantics.
JNIEXPORT jobject JNICALL
Java_org_kensmath_gopack_GOPackNative_computeRandomDisc(
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
Java_org_kensmath_gopack_GOPackNative_nativeVersion(JNIEnv* env, jclass /*clazz*/) {
    return env->NewStringUTF(
        "gopack-cpp 0.1.0 (max-pack mode ported; polygonal/rectangle modes pending)");
}

} // extern "C"
