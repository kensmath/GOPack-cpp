// GOPack C++ -- standalone command-line front end.
//
// Usage: gopack <input.p> -o <output.p> [--passes 20] [--eucl-out]
//        gopack <input.p> -o <output.p> --polygon [--corners v1,v2,v3,v4] [--angles a1,a2,a3,a4]
//
// Loads a *.p triangulation/packing file, runs the maximal-packing (or, with
// --polygon, polygonal/rectangle) riffle loop, and writes the resulting
// packing back out. This is a thin wrapper around gopack::Packer; it exists
// so the core engine can be exercised as a standalone Windows/macOS
// executable, e.g. for scripting or for a Java caller that prefers
// subprocess invocation over JNI (see jni/cpp for the in-process JNI
// alternative, which wraps the same Packer class).
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#if defined(_MSC_VER) && defined(_DEBUG)
#include <crtdbg.h>
#endif

#include "gopack/Packer.h"

namespace {

void printUsage(const char* argv0) {
    std::fprintf(stderr,
                 "Usage: %s <input.p> -o <output.p> [--passes N] [--eucl-out]\n"
                 "       %s <input.p> -o <output.p> --polygon [--corners v1,v2,v3,...] "
                 "[--angles a1,a2,a3,...]\n"
                 "  --polygon           use polygonal/rectangle packing mode (setMode mode 2)\n"
                 "                      instead of the default maximal-packing mode (mode 1)\n"
                 "  --corners v1,v2,... comma-separated 1-indexed boundary vertices to use as\n"
                 "                      polygon corners, in counterclockwise order; if omitted,\n"
                 "                      corners are chosen automatically (see setMode.m)\n"
                 "  --angles a1,a2,...  comma-separated target corner angles (radians), matching\n"
                 "                      --corners in count and order; if omitted, corners get\n"
                 "                      equal angles (e.g. pi/2 each for a 4-corner rectangle)\n",
                 argv0, argv0);
}

// Parses a comma-separated list of integers, e.g. "3,17,42".
std::vector<gopack::Index> parseIndexList(const std::string& s) {
    std::vector<gopack::Index> out;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(std::stoi(tok));
    }
    return out;
}

// Parses a comma-separated list of doubles, e.g. "1.5708,1.5708,1.5708,1.5708".
std::vector<gopack::Scalar> parseScalarList(const std::string& s) {
    std::vector<gopack::Scalar> out;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) out.push_back(std::stod(tok));
    }
    return out;
}

} // namespace

int main(int argc, char** argv) {
#if defined(_MSC_VER) && defined(_DEBUG)
    // By default, an MSVC Debug build's checked-iterator failures ("vector
    // subscript out of range", "vector iterator not dereferenceable", etc.)
    // pop up a blocking modal assertion dialog box rather than printing to
    // the console. That dialog has no text you can copy from a terminal --
    // which is exactly what happened on a previous run of this program: the
    // message was seen but couldn't be captured. Redirect all CRT
    // assert/error reporting to stderr instead, so the exact expression,
    // file, and line print to the console (and can be piped to a log file)
    // instead of blocking on an invisible-to-the-terminal dialog.
    _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDERR);
    _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
    _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDERR);
#endif
    if (argc < 2) {
        printUsage(argv[0]);
        return 2;
    }

    std::string inputPath = argv[1];
    std::string outputPath;
    int passes = 20;
    bool euclOut = false;
    bool polygonMode = false;
    std::vector<gopack::Index> corners;
    std::vector<gopack::Scalar> angles;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-o" && i + 1 < argc) {
            outputPath = argv[++i];
        } else if (arg == "--passes" && i + 1 < argc) {
            passes = std::stoi(argv[++i]);
        } else if (arg == "--eucl-out") {
            euclOut = true;
        } else if (arg == "--polygon") {
            polygonMode = true;
        } else if (arg == "--corners" && i + 1 < argc) {
            corners = parseIndexList(argv[++i]);
        } else if (arg == "--angles" && i + 1 < argc) {
            angles = parseScalarList(argv[++i]);
        } else {
            std::fprintf(stderr, "Unrecognized argument: %s\n", arg.c_str());
            printUsage(argv[0]);
            return 2;
        }
    }

    try {
        gopack::Packer packer;
        if (packer.readpack(inputPath) <= 0) {
            std::fprintf(stderr, "gopack: failed to read %s\n", inputPath.c_str());
            return 1;
        }

        if (polygonMode) {
            // Note: no need to call indxMatrices() again here -- mode only
            // changes boundary layout (vAims/corners/sides), not which
            // vertices are interior, so the layoutVerts/rimVerts/tranI-etc.
            // data indxMatrices() built inside readpack() is still valid.
            if (packer.setMode(2, corners, angles) < 0) {
                std::fprintf(stderr, "gopack: failed to set polygonal packing mode\n");
                return 1;
            }
        } else {
            if (packer.setMode(1) < 0) {
                std::fprintf(stderr, "gopack: failed to set max-pack mode\n");
                return 1;
            }
        }

        gopack::RiffleResult result = packer.riffle(passes);
        if (result.cycles < 0) {
            std::fprintf(stderr, "gopack: riffle failed\n");
            return 1;
        }

        packer.packStatus();
        if (polygonMode && packer.corners.size() == 4) {
            std::fprintf(stdout, "  Aspect ratio: %f\n", packer.getAspect());
        }

        if (!outputPath.empty()) {
            if (euclOut) {
                packer.writeEucl(outputPath);
            } else {
                packer.writepack(outputPath);
            }
        }

        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "gopack: %s\n", e.what());
        return 1;
    }
}
