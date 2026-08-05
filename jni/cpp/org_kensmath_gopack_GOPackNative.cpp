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

JNIEXPORT jstring JNICALL
Java_org_kensmath_gopack_GOPackNative_nativeVersion(JNIEnv* env, jclass /*clazz*/) {
    return env->NewStringUTF(
        "gopack-cpp 0.1.0 (max-pack mode ported; polygonal/rectangle modes pending)");
}

} // extern "C"
