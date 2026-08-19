// GOPack C++ -- standalone command-line front end.
//
// Usage: gopack <input.p> -o <output.p> [--passes 200] [--eucl-out]
//        gopack <input.p> -o <output.p> --polygon [--corners v1,v2,v3,v4] [--angles a1,a2,a3,a4]
//        gopack --random-disc N -o <output.p> [--passes 200] [--eucl-out]
//        gopack --random-sphere N -o <output.p> [--passes 200]
//        gopack --random-square N -o <output.p> [--passes 200] [--eucl-out]
//        gopack --random-rectangle N[,aspect[,bdryN]] -o <output.p> [--passes 200] [--eucl-out]
//        gopack --random-tri intN,bdryN --graph x1,y1,x2,y2,... -o <output.p> [--cent cx,cy]
//               [--passes 200] [--eucl-out]
//
// Loads a *.p triangulation/packing file, or -- with one of the
// --random-disc/--random-sphere/--random-square/--random-rectangle/
// --random-tri flags -- generates a fresh "geometrically random"
// triangulation instead of reading one (see gopack/RandomGen.h and
// Packer::randomDisc()/randomSphere()/randomRectangle()/randomSquare()/
// randomTri(), the C++ port of randomDisc.m/randomSphere.m/
// randomRectangle.m/randomSquare.m/randomTri.m -- only available when this
// binary was built with GOPACK_BUILD_RANDOM_GEN, the default). --random-tri
// is the general case: it fills an arbitrary user-supplied closed polygonal
// region ('--graph', a flat x,y coordinate list) rather than a fixed disc/
// square/rectangle shape -- see Packer::randomTri()'s doc comment in
// Packer.h for the full semantics (in particular: unlike --random-rectangle,
// this leaves mode at 1/max-pack, not polygonal mode, so --polygon/
// --corners/--angles still apply afterward if you want polygonal mode on
// the result). Either way, this then runs the maximal-packing (or, with
// --polygon, polygonal/rectangle) riffle loop, and writes the resulting
// packing back out. This is a thin wrapper around gopack::Packer; it exists
// so the core engine can be exercised as a standalone Windows/macOS
// executable, e.g. for scripting or for a Java caller that prefers
// subprocess invocation over JNI (see jni/cpp for the in-process JNI
// alternative, which wraps the same Packer class and also exposes all five
// --random-* generators as of the computeRandom*() bridge methods).
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
#ifdef GOPACK_HAVE_RANDOM_GEN
                 "       %s --random-disc N -o <output.p> [--passes N] [--eucl-out]\n"
                 "       %s --random-sphere N -o <output.p> [--passes N]\n"
                 "       %s --random-square N -o <output.p> [--passes N] [--eucl-out]\n"
                 "       %s --random-rectangle N[,aspect[,bdryN]] -o <output.p> [--passes N] "
                 "[--eucl-out]\n"
                 "       %s --random-tri intN,bdryN --graph x1,y1,x2,y2,... -o <output.p> "
                 "[--cent cx,cy] [--passes N] [--eucl-out]\n"
#endif
                 "  --polygon           use polygonal/rectangle packing mode (setMode mode 2)\n"
                 "                      instead of the default maximal-packing mode (mode 1)\n"
                 "  --corners v1,v2,... comma-separated 1-indexed boundary vertices to use as\n"
                 "                      polygon corners, in counterclockwise order; if omitted,\n"
                 "                      corners are chosen automatically (see setMode.m)\n"
                 "  --angles a1,a2,...  comma-separated target corner angles (radians), matching\n"
                 "                      --corners in count and order; if omitted, corners get\n"
                 "                      equal angles (e.g. pi/2 each for a 4-corner rectangle)\n"
#ifdef GOPACK_HAVE_RANDOM_GEN
                 "  --random-disc N            generate a random triangulation of the unit disc\n"
                 "                             with N total points instead of reading <input.p>\n"
                 "  --random-sphere N          generate a random triangulation of the sphere\n"
                 "                             with N points instead of reading <input.p>\n"
                 "  --random-square N          generate a random triangulation of the unit\n"
                 "                             square with N total points (already in polygonal\n"
                 "                             mode; --polygon/--corners/--angles are ignored)\n"
                 "  --random-rectangle N[,aspect[,bdryN]]\n"
                 "                             generate a random triangulation of the rectangle\n"
                 "                             [-aspect,aspect]x[-1,1] (default aspect 1) with N\n"
                 "                             interior points (already in polygonal mode;\n"
                 "                             --polygon/--corners/--angles are ignored)\n"
                 "  --random-tri intN,bdryN    generate a random triangulation of an arbitrary\n"
                 "                             closed polygonal region (see --graph) with intN\n"
                 "                             interior and bdryN boundary points; stays in\n"
                 "                             max-pack mode (--polygon/--corners/--angles still\n"
                 "                             apply afterward if you want polygonal mode)\n"
                 "  --graph x1,y1,x2,y2,...    required with --random-tri: the closed boundary\n"
                 "                             polygon's vertices, as a flat x,y coordinate list\n"
                 "                             (do not repeat the first point at the end)\n"
                 "  --cent cx,cy               optional with --random-tri: a point inside --graph\n"
                 "                             to use as the packing's alpha vertex; if omitted,\n"
                 "                             or the point isn't inside --graph, alpha is chosen\n"
                 "                             automatically\n"
#endif
                 ,
#ifdef GOPACK_HAVE_RANDOM_GEN
                 argv0, argv0, argv0, argv0, argv0, argv0, argv0);
#else
                 argv0, argv0);
#endif
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

#ifdef GOPACK_HAVE_RANDOM_GEN
// Parses a flat "x1,y1,x2,y2,..." coordinate list (as used by --graph/
// --cent) into gopack::Complex points. Returns false (leaving 'out'
// unspecified) if the count of numbers is odd -- every point needs both an
// x and a y.
bool parseComplexList(const std::string& s, std::vector<gopack::Complex>& out) {
    std::vector<gopack::Scalar> flat = parseScalarList(s);
    if (flat.size() % 2 != 0) return false;
    out.clear();
    out.reserve(flat.size() / 2);
    for (size_t i = 0; i + 1 < flat.size(); i += 2) {
        out.emplace_back(flat[i], flat[i + 1]);
    }
    return true;
}
#endif

// What the packing comes from: an existing *.p file, or one of the
// --random-* generators (only reachable when GOPACK_HAVE_RANDOM_GEN is
// defined -- see the arg1 dispatch in main()).
enum class Source { File, RandomDisc, RandomSphere, RandomSquare, RandomRectangle, RandomTri };

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

    Source source = Source::File;
    std::string inputPath;
    gopack::Index randomN = 0;
    gopack::Scalar randomAspect = 1.0;
    gopack::Index randomBdryN = -1;
#ifdef GOPACK_HAVE_RANDOM_GEN
    std::vector<gopack::Complex> randomGraph;
    bool haveCent = false;
    gopack::Complex randomCent;
#endif
    int flagsStart = 2;

    const std::string arg1 = argv[1];
#ifdef GOPACK_HAVE_RANDOM_GEN
    if (arg1 == "--random-disc" || arg1 == "--random-sphere" || arg1 == "--random-square") {
        if (argc < 3) {
            std::fprintf(stderr, "%s requires an integer point count\n", arg1.c_str());
            printUsage(argv[0]);
            return 2;
        }
        randomN = std::stoi(argv[2]);
        source = (arg1 == "--random-disc")
                     ? Source::RandomDisc
                     : (arg1 == "--random-sphere") ? Source::RandomSphere : Source::RandomSquare;
        flagsStart = 3;
    } else if (arg1 == "--random-rectangle") {
        if (argc < 3) {
            std::fprintf(stderr, "--random-rectangle requires N[,aspect[,bdryN]]\n");
            printUsage(argv[0]);
            return 2;
        }
        std::vector<std::string> parts;
        std::istringstream ss(argv[2]);
        std::string tok;
        while (std::getline(ss, tok, ',')) parts.push_back(tok);
        if (parts.empty() || parts[0].empty()) {
            std::fprintf(stderr, "--random-rectangle requires N[,aspect[,bdryN]]\n");
            printUsage(argv[0]);
            return 2;
        }
        randomN = std::stoi(parts[0]);
        if (parts.size() > 1 && !parts[1].empty()) randomAspect = std::stod(parts[1]);
        if (parts.size() > 2 && !parts[2].empty()) randomBdryN = std::stoi(parts[2]);
        source = Source::RandomRectangle;
        flagsStart = 3;
    } else if (arg1 == "--random-tri") {
        if (argc < 3) {
            std::fprintf(stderr, "--random-tri requires intN,bdryN\n");
            printUsage(argv[0]);
            return 2;
        }
        std::vector<std::string> parts;
        std::istringstream ss(argv[2]);
        std::string tok;
        while (std::getline(ss, tok, ',')) parts.push_back(tok);
        if (parts.size() != 2 || parts[0].empty() || parts[1].empty()) {
            std::fprintf(stderr, "--random-tri requires intN,bdryN (both required)\n");
            printUsage(argv[0]);
            return 2;
        }
        randomN = std::stoi(parts[0]);
        randomBdryN = std::stoi(parts[1]);
        source = Source::RandomTri;
        flagsStart = 3;
    } else
#endif
    {
        inputPath = arg1;
        flagsStart = 2;
    }

    std::string outputPath;
    // Bug-fix-adjacent tuning, not a correctness issue: 20 was this CLI's
    // original default (an arbitrary choice made when this file was first
    // written, not inherited from anywhere in MATLAB), and it turned out
    // too low to fully converge some small/irregular "geometrically
    // random" triangulations (e.g. --random-disc with N in the tens of
    // points) -- riffle() stopping before visErr reaches continueRiffle()'s
    // own cutval=0.01 convergence threshold leaves neighboring circles'
    // radii/centers not yet satisfying the tangency constraints, which is
    // visible as overlapping circles once loaded in CirclePack. Bumping
    // this default costs nothing for inputs that already converge quickly
    // (a real ~25-30k vertex packing, for instance): continueRiffle()'s own
    // loop already exits the moment maxVis <= cutval, well before passNum
    // passes if the packing gets there sooner -- this only helps inputs
    // that were genuinely stopping short before.
    int passes = 200;
    bool euclOut = false;
    bool polygonMode = false;
    std::vector<gopack::Index> corners;
    std::vector<gopack::Scalar> angles;

    for (int i = flagsStart; i < argc; ++i) {
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
#ifdef GOPACK_HAVE_RANDOM_GEN
        } else if (arg == "--graph" && i + 1 < argc) {
            if (!parseComplexList(argv[++i], randomGraph)) {
                std::fprintf(stderr,
                              "--graph requires an even count of numbers (x1,y1,x2,y2,...)\n");
                return 2;
            }
        } else if (arg == "--cent" && i + 1 < argc) {
            std::vector<gopack::Complex> centPts;
            if (!parseComplexList(argv[++i], centPts) || centPts.size() != 1) {
                std::fprintf(stderr, "--cent requires exactly one point (cx,cy)\n");
                return 2;
            }
            randomCent = centPts[0];
            haveCent = true;
#endif
        } else {
            std::fprintf(stderr, "Unrecognized argument: %s\n", arg.c_str());
            printUsage(argv[0]);
            return 2;
        }
    }

#ifdef GOPACK_HAVE_RANDOM_GEN
    if (source == Source::RandomTri && randomGraph.size() < 3) {
        std::fprintf(stderr, "--random-tri requires --graph with at least 3 points\n");
        printUsage(argv[0]);
        return 2;
    }
#endif

    try {
        gopack::Packer packer;
        // Random-rectangle/random-square already leave the packer in
        // polygonal mode (setMode(2, ...) with corners chosen from the
        // generated rectangle -- see Packer::randomRectangle()), so the
        // usual setMode() call below (and --polygon/--corners/--angles) is
        // skipped for those two sources.
        bool modeAlreadySet = false;

#ifdef GOPACK_HAVE_RANDOM_GEN
        switch (source) {
            case Source::RandomDisc:
                packer = gopack::Packer::randomDisc(randomN);
                break;
            case Source::RandomSphere:
                packer = gopack::Packer::randomSphere(randomN);
                break;
            case Source::RandomSquare:
                packer = gopack::Packer::randomSquare(randomN);
                modeAlreadySet = true;
                break;
            case Source::RandomRectangle:
                packer = gopack::Packer::randomRectangle(randomN, randomAspect, randomBdryN);
                modeAlreadySet = true;
                break;
            case Source::RandomTri:
                // Leaves mode at 1 (max-pack), not polygonal mode -- unlike
                // RandomSquare/RandomRectangle, this is a generic region with
                // no "corners" concept, so modeAlreadySet stays false and the
                // usual setMode()/--polygon handling below still applies.
                packer = gopack::Packer::randomTri(randomN, randomBdryN, randomGraph,
                                                    haveCent ? &randomCent : nullptr);
                break;
            case Source::File:
                if (packer.readpack(inputPath) <= 0) {
                    std::fprintf(stderr, "gopack: failed to read %s\n", inputPath.c_str());
                    return 1;
                }
                break;
        }
#else
        (void)source;
        (void)randomN;
        (void)randomAspect;
        (void)randomBdryN;
        if (packer.readpack(inputPath) <= 0) {
            std::fprintf(stderr, "gopack: failed to read %s\n", inputPath.c_str());
            return 1;
        }
#endif

        if (packer.nodeCount <= 0) {
            std::fprintf(stderr, "gopack: failed to generate a usable packing\n");
            return 1;
        }

        if (!modeAlreadySet) {
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
        }

        gopack::RiffleResult result = packer.riffle(passes);
        if (result.cycles < 0) {
            std::fprintf(stderr, "gopack: riffle failed\n");
            return 1;
        }

        packer.packStatus();
        if (packer.mode == 2 && packer.corners.size() == 4) {
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
