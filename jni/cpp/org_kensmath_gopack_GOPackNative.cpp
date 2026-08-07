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

        gopack::RiffleResult result = packer.riffle(maxPasses > 0 ? maxPasses : 20);
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

        gopack::RiffleResult result = packer.riffle(maxPasses > 0 ? maxPasses : 20);
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

JNIEXPORT jstring JNICALL
Java_org_kensmath_gopack_GOPackNative_nativeVersion(JNIEnv* env, jclass /*clazz*/) {
    return env->NewStringUTF(
        "gopack-cpp 0.1.0 (max-pack mode ported; polygonal/rectangle modes pending)");
}

} // extern "C"
