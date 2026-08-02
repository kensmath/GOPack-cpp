# JNI shared library target.
#
# Requires JAVA_HOME (or -DJNI_INCLUDE_DIRS=... / -DJava_INCLUDE_DIRS=...) to
# be set so find_package(JNI) can locate <jni.h>. This target is optional
# (guarded by GOPACK_BUILD_JNI in the top-level CMakeLists.txt) precisely so
# the CLI/core can still be built on machines without a JDK installed.
find_package(JNI)

if(NOT JNI_FOUND)
    message(WARNING
        "GOPACK_BUILD_JNI is ON but JNI headers were not found (JAVA_HOME not "
        "set?). Skipping the JNI shared library target.")
    return()
endif()

add_library(gopack_jni SHARED
    org_kensmath_gopack_GOPackNative.cpp
)

target_include_directories(gopack_jni PRIVATE ${JNI_INCLUDE_DIRS})
target_link_libraries(gopack_jni PRIVATE gopack_core ${JNI_LIBRARIES})

# Produces libgopack_jni.dylib on macOS, gopack_jni.dll on Windows (MSVC/MinGW
# both honor this via the SHARED target + platform default naming), and
# libgopack_jni.so on Linux.
set_target_properties(gopack_jni PROPERTIES
    OUTPUT_NAME "gopack_jni"
)

install(TARGETS gopack_jni
    RUNTIME DESTINATION bin      # .dll lands here on Windows
    LIBRARY DESTINATION lib      # .dylib / .so land here
)
