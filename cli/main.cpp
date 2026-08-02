// GOPack C++ -- standalone command-line front end.
//
// Usage: gopack <input.p> -o <output.p> [--passes 20] [--eucl-out]
//
// Loads a *.p triangulation/packing file, runs the maximal-packing riffle
// loop, and writes the resulting packing back out. This is a thin wrapper
// around gopack::Packer; it exists so the core engine can be exercised as a
// standalone Windows/macOS executable, e.g. for scripting or for a Java
// caller that prefers subprocess invocation over JNI (see jni/cpp for the
// in-process JNI alternative, which wraps the same Packer class).
#include <cstdio>
#include <string>

#include "gopack/Packer.h"

namespace {

void printUsage(const char* argv0) {
    std::fprintf(stderr,
                 "Usage: %s <input.p> -o <output.p> [--passes N] [--eucl-out]\n",
                 argv0);
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 2;
    }

    std::string inputPath = argv[1];
    std::string outputPath;
    int passes = 20;
    bool euclOut = false;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-o" && i + 1 < argc) {
            outputPath = argv[++i];
        } else if (arg == "--passes" && i + 1 < argc) {
            passes = std::stoi(argv[++i]);
        } else if (arg == "--eucl-out") {
            euclOut = true;
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

        if (packer.setMode(1) < 0) {
            std::fprintf(stderr, "gopack: failed to set max-pack mode\n");
            return 1;
        }

        gopack::RiffleResult result = packer.riffle(passes);
        if (result.cycles < 0) {
            std::fprintf(stderr, "gopack: riffle failed\n");
            return 1;
        }

        packer.packStatus();

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
