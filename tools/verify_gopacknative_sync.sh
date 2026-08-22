#!/usr/bin/env bash
# Verifies that jni/java/JNI/GOPackNative.java's declared native methods
# exactly match jni/cpp/JNI_GOPackNative.cpp's actual exported symbols --
# both which methods exist, and their JNI return-type category (the
# object-array vs. flat-array distinction that caused the
# EXCEPTION_ACCESS_VIOLATION crash once, and the double[]/double[][]
# regression HANDOFFmaxpackneedscenters.md reported).
#
# Run this BEFORE accepting any regenerated/resynced copy of
# GOPackNative.java into either repo (GOPack-cpp's own copy, or
# CirclePack's vendored copy) -- a mismatch here means the Java side and
# the native side have silently drifled apart, which links and runs
# without any compiler/linker error and only fails at runtime (or, for a
# genuinely missing method, fails to compile the CALLER, not this file).
#
# Requires only `javac` and `grep`/`sed` -- no JDK native headers, no C++
# compiler, no Eigen -- so it runs anywhere either repo is checked out.
#
# Usage: verify_gopacknative_sync.sh <path-to-GOPack-cpp-repo-root>
set -euo pipefail

REPO="${1:?usage: $0 <path-to-GOPack-cpp-repo-root>}"
JAVA_SRC_DIR="$REPO/jni/java/JNI"
CPP_FILE="$REPO/jni/cpp/JNI_GOPackNative.cpp"

for f in "$JAVA_SRC_DIR/GOPackNative.java" "$JAVA_SRC_DIR/RandomComplexResult.java" \
         "$JAVA_SRC_DIR/GOPackException.java" "$CPP_FILE"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: expected file not found: $f" >&2
        exit 2
    fi
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/JNI" "$WORK/out"
cp "$JAVA_SRC_DIR"/*.java "$WORK/JNI/"

echo "== Compiling Java side and generating JNI header =="
( cd "$WORK" && javac -h . -d out JNI/*.java )

JAVA_HEADER="$WORK/JNI_GOPackNative.h"
if [ ! -f "$JAVA_HEADER" ]; then
    echo "ERROR: javac -h produced no JNI_GOPackNative.h -- does GOPackNative.java" >&2
    echo "       declare any native methods at all?" >&2
    exit 2
fi

# Extract "name -> return-type-and-arg-descriptor" from the javac-generated
# header (ground truth for what the JAVA side declares).
extract_java_sigs() {
    awk '
        /^JNIEXPORT/ { rettype = $2 }
        /^  \(JNIEnv/ {
            name = prevname
            print name " " rettype " " $0
        }
        /^JNIEXPORT.*Java_JNI_GOPackNative_/ {
            match($0, /Java_JNI_GOPackNative_[A-Za-z0-9_]+/)
            prevname = substr($0, RSTART, RLENGTH)
        }
    ' "$JAVA_HEADER" | sort
}

# Extract "name -> return-type" from the actual .cpp implementations
# (ground truth for what the NATIVE side actually returns).
extract_cpp_sigs() {
    grep -oP '^(JNIEXPORT\s+\S+\s+JNICALL\s*\n?)?Java_JNI_GOPackNative_[A-Za-z0-9_]+' "$CPP_FILE" >/dev/null 2>&1 || true
    # More robust: scan JNIEXPORT lines, which are always immediately
    # followed (possibly on the same or next line) by the Java_... symbol.
    awk '
        /^JNIEXPORT/ {
            rettype = $2
            getline nextline
            if (match(nextline, /Java_JNI_GOPackNative_[A-Za-z0-9_]+/)) {
                name = substr(nextline, RSTART, RLENGTH)
                print name " " rettype
            }
        }
    ' "$CPP_FILE" | sort
}

echo "== Cross-checking method names and return-type categories =="
JAVA_NAMES="$(extract_java_sigs | awk '{print $1}' | sort -u)"
CPP_NAMES="$(extract_cpp_sigs | awk '{print $1}' | sort -u)"

FAIL=0

MISSING_IN_JAVA="$(comm -23 <(echo "$CPP_NAMES") <(echo "$JAVA_NAMES"))"
if [ -n "$MISSING_IN_JAVA" ]; then
    echo "FAIL: .cpp exports these but GOPackNative.java declares no matching native method:"
    echo "$MISSING_IN_JAVA" | sed 's/^/  /'
    FAIL=1
fi

EXTRA_IN_JAVA="$(comm -13 <(echo "$CPP_NAMES") <(echo "$JAVA_NAMES"))"
if [ -n "$EXTRA_IN_JAVA" ]; then
    echo "FAIL: GOPackNative.java declares these native methods but .cpp exports no matching symbol:"
    echo "$EXTRA_IN_JAVA" | sed 's/^/  /'
    FAIL=1
fi

echo
echo "== Cross-checking return-type category per method (object-array vs. flat-array is the dangerous mismatch) =="
while read -r name javaret _; do
    cppret="$(extract_cpp_sigs | awk -v n="$name" '$1==n {print $2}')"
    [ -z "$cppret" ] && continue  # already reported as missing above
    # Normalize to a coarse category: jobjectArray/jobject vs. everything
    # else (jdoubleArray/jintArray/jstring/jdouble/jint/...) -- this is
    # exactly the axis that caused the original crash and this regression.
    javacat="$javaret"; cppcat="$cppret"
    case "$javaret" in jobjectArray|jobject) javacat="OBJECT";; *) javacat="PRIMITIVE";; esac
    case "$cppret" in jobjectArray|jobject) cppcat="OBJECT";; *) cppcat="PRIMITIVE";; esac
    if [ "$javacat" != "$cppcat" ]; then
        echo "FAIL: $name -- Java header says $javaret ($javacat), .cpp returns $cppret ($cppcat)"
        FAIL=1
    fi
done <<< "$(extract_java_sigs)"

if [ "$FAIL" -eq 0 ]; then
    echo "OK: all native methods match between GOPackNative.java and JNI_GOPackNative.cpp."
    exit 0
else
    echo
    echo "GOPackNative.java is NOT in sync with JNI_GOPackNative.cpp -- do not commit/accept" >&2
    echo "this copy until the mismatches above are resolved." >&2
    exit 1
fi
