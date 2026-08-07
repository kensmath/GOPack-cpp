#include "gopack/Packer.h"

#include <cctype>
#include <cstdio>
#include <fstream>
#include <istream>
#include <sstream>

#include "gopack/Geometry.h"

namespace gopack {

namespace {

// trimFilename.m -- base name without directory or extension.
std::string trimFilename(const std::string& path) {
    size_t slash = path.find_last_of("/\\");
    std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    size_t dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    return base;
}

bool startsWith(const std::string& s, const char* prefix) {
    size_t len = std::char_traits<char>::length(prefix);
    return s.size() >= len && s.compare(0, len, prefix) == 0;
}

// Attempts to read one integer; on failure, rewinds the stream to just
// before the attempt (clearing any fail state) so the un-consumed token is
// left for the caller to read as something else. This mirrors MATLAB's
// fscanf(fid,'%d',N), which reads as many of the N requested integers as it
// can and simply stops -- without erroring or losing its place -- the
// moment it hits a non-numeric token. A plain `in >> x` does not have that
// "stop gracefully" behavior, which matters here: some *.p files write
// "ALPHA/GAMMA: a g" (2 values) while others write
// "ALPHA/BETA/GAMMA: a b g" (3 values), and the reader has to accept both.
bool tryReadInt(std::istream& in, Index& out) {
    auto pos = in.tellg();
    if (in >> out) return true;
    in.clear();
    in.seekg(pos);
    return false;
}

} // namespace

Index Packer::readpack(const std::string& fname) {
    // Open in binary mode deliberately. On Windows, a text-mode ifstream
    // translates "\r\n" -> "\n" on read, and tellg()/seekg() on a text-mode
    // stream return/accept opaque positions that are only reliably valid at
    // very specific points -- they are NOT guaranteed to be plain byte
    // offsets, and round-tripping them across a formatted extraction
    // attempt (as tryReadInt does below, to support both the 2- and
    // 3-value ALPHA/GAMMA and ALPHA/BETA/GAMMA variants) can silently land
    // the read position at the wrong place. That desyncs every token read
    // after the rewind -- in practice this was observed to skip the
    // "BOUQUET:"/"FLOWERS:" keyword entirely, leaving `flowers` empty and
    // crashing complexCount() on its very first access. Binary mode makes
    // tellg()/seekg() plain, reliable byte offsets on every platform; the
    // token/getline-based parsing below already treats '\r' as whitespace
    // (std::isspace includes it), so a stray trailing '\r' per line from
    // CRLF files is harmless.
    std::ifstream in(fname, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "Failed to read packing file %s\n", fname.c_str());
        return 0;
    }

    int nodeCase = 1;
    int geomFlag = 0; // -1 hyp, 0 eucl, 1 sph
    std::vector<Scalar> newRadii;
    std::vector<Complex> newCenters;
    bool haveNewRadii = false, haveNewCenters = false;

    std::string beginStr;
    in >> beginStr;

    if (startsWith(beginStr, "CHECKCOUNT:")) {
        nodeCase = 0;
        Index ckcount = 0;
        in >> ckcount;
        if (ckcount != nodeCount) {
            std::fprintf(stderr,
                         "Error: this is a CHECKCOUNT type file, count does not match\n");
            return 0;
        }
    } else if (startsWith(beginStr, "NODECOUNT:")) {
        cleanse();
        in >> nodeCount;
        fileName = trimFilename(fname);
    } else if (startsWith(beginStr, "OFF")) {
        std::fprintf(stderr, "TODO: readpack in OFF format; code not yet done.\n");
        return 0;
    } else {
        // Attempt to read as a bare triangle list (v1 v2 v3 per face).
        // parse_triangles.m has not been ported yet; report and bail out
        // rather than guess, matching this file's other "not yet done" paths.
        std::fprintf(stderr,
                     "Packer::readpack: input does not start with NODECOUNT:/CHECKCOUNT:, "
                     "and the bare-triangle-list reader (parse_triangles.m) has not been "
                     "ported yet.\n");
        return 0;
    }

    bool done = false;
    while (!done) {
        std::string tok;
        if (!(in >> tok)) break;

        if (startsWith(tok, "ALPH")) {
            // Two file variants exist: "ALPHA/BETA/GAMMA: a b g" (3 values,
            // beta unused) and "ALPHA/GAMMA: a g" (2 values). Read as many
            // as are actually present, matching the MATLAB fscanf behavior
            // described above tryReadInt's definition.
            Index a = 0, b = 0, c = 0;
            bool haveA = tryReadInt(in, a);
            bool haveB = haveA && tryReadInt(in, b);
            bool haveC = haveB && tryReadInt(in, c);
            if (haveA) alpha = a;
            if (haveC) {
                gamma = c; // 3 values: alpha, beta (unused), gamma
            } else if (haveB) {
                gamma = b; // 2 values: alpha, gamma
            }
        } else if (startsWith(tok, "GEOM")) { // GEOMETRY:
            std::string g;
            in >> g;
            if (g.find("hyp") != std::string::npos) {
                geomFlag = -1;
            } else if (g.find("sph") != std::string::npos) {
                geomFlag = 1;
            } else {
                geomFlag = 0;
            }
            if (nodeCase == 1) hes = static_cast<Geometry>(geomFlag);
        } else if ((startsWith(tok, "FLOW") || startsWith(tok, "BOUQ")) && nodeCase == 1) {
            std::string restOfLine;
            std::getline(in, restOfLine); // consume rest of "FLOWERS:" line

            std::vector<std::vector<Index>> parsedFlowers(static_cast<size_t>(nodeCount) + 1);
            for (Index j = 1; j <= nodeCount; ++j) {
                std::string line;
                while (std::getline(in, line)) {
                    // skip blank/whitespace-only lines
                    bool blank = true;
                    for (char c : line) {
                        if (!std::isspace(static_cast<unsigned char>(c))) { blank = false; break; }
                    }
                    if (!blank) break;
                }
                std::istringstream ls(line);
                std::vector<Scalar> nums;
                Scalar val;
                while (ls >> val) nums.push_back(val);

                if (nums.size() < 3) {
                    std::fprintf(stderr,
                                 "Packer::readpack: malformed FLOWERS line for vertex %d\n", j);
                    return 0;
                }
                // nums[1] is the file's own declared petal count; from here
                // on only the flower list itself (nums[2..end]) is kept, and
                // ingestFlowers() derives vNum as flower.size()-1. For any
                // well-formed *.p file the two always agree (that's what
                // "well-formed" means here), so this isn't a behavior change
                // in practice -- it just avoids readpack() and the new
                // in-memory loadComplex() entry point needing two different
                // notions of vNum.
                std::vector<Index> flower;
                flower.reserve(nums.size() - 2);
                for (size_t idx = 2; idx < nums.size(); ++idx) {
                    flower.push_back(static_cast<Index>(nums[idx]));
                }
                parsedFlowers[static_cast<size_t>(j)] = std::move(flower);
            }

            if (ingestFlowers(nodeCount, parsedFlowers) == 0) {
                return 0;
            }
        } else if (startsWith(tok, "RADI")) { // RADII:
            newRadii.assign(static_cast<size_t>(nodeCount) + 1, 0.0);
            for (Index v = 1; v <= nodeCount; ++v) in >> newRadii[v];
            haveNewRadii = true;
        } else if (startsWith(tok, "CENT")) { // CENTERS:
            newCenters.assign(static_cast<size_t>(nodeCount) + 1, Complex(0.0, 0.0));
            for (Index v = 1; v <= nodeCount; ++v) {
                Scalar x = 0, y = 0;
                in >> x >> y;
                newCenters[v] = Complex(x, y);
            }
            haveNewCenters = true;
        } else if (startsWith(tok, "VERT")) { // VERT_LIST
            vlist.clear();
            Scalar v;
            for (Index j = 0; j < nodeCount && (in >> v); ++j) {
                vlist.push_back(static_cast<Index>(v));
            }
        } else if (startsWith(tok, "END")) {
            done = true;
        }
    }

    if (nodeCase == 1) {
        // Fresh load (as opposed to a CHECKCOUNT: update to an
        // already-loaded packing, handled in the else branch below): this is
        // exactly what the new loadComplex() entry point also needs to do
        // once it has a complex in hand, so it's factored out into
        // finalizeComplex() and shared rather than duplicated.
        const std::vector<Scalar>* radiiPtr = haveNewRadii ? &newRadii : nullptr;
        const std::vector<Complex>* centersPtr = haveNewCenters ? &newCenters : nullptr;
        if (finalizeComplex(radiiPtr, centersPtr, nullptr, fileName) == 0) {
            return 0;
        }
    } else {
        // CHECKCOUNT-type file: update radii/centers on the already-loaded
        // packing in place, without touching combinatorics/alpha/vAims (that
        // state belongs to the packing that was already loaded by an
        // earlier readpack()/loadComplex() call).
        if (geomFlag == 0 && haveNewRadii) radii = newRadii;
        if (geomFlag == 0 && haveNewCenters) centers = newCenters;
        if (geomFlag < 0 && haveNewRadii && haveNewCenters) {
            for (Index v = 1; v <= nodeCount; ++v) {
                auto [ez, er] = geom::hToEData(newCenters[v], newRadii[v]);
                radii[v] = er;
                centers[v] = ez;
            }
        }
        localcenters = centers;
        localradii = radii;
        mode = 1;
        indxMatrices();
    }

    return nodeCount;
}

Index Packer::ingestFlowers(Index nodeCountIn, const std::vector<std::vector<Index>>& flowersIn) {
    nodeCount = nodeCountIn;
    flowers.assign(static_cast<size_t>(nodeCount) + 1, {});
    vNum.assign(static_cast<size_t>(nodeCount) + 1, 0);
    bdryFlags.assign(static_cast<size_t>(nodeCount) + 1, 0);
    bdryCount = 0;

    for (Index j = 1; j <= nodeCount; ++j) {
        const std::vector<Index>& flower = flowersIn[static_cast<size_t>(j)];
        if (flower.size() < 2) {
            std::fprintf(stderr,
                         "Packer::ingestFlowers: flower for vertex %d has fewer than 2 "
                         "entries (need at least a single petal, closed or open)\n", j);
            return 0;
        }
        vNum[j] = static_cast<Index>(flower.size()) - 1;
        if (flower.front() != flower.back()) { // bdry vertex
            bdryFlags[j] = 1;
            bdryCount++;
        }
        flowers[j] = flower;
    }

    intCount = nodeCount - bdryCount;
    Index totvNum = 0;
    for (Index v = 1; v <= nodeCount; ++v) totvNum += vNum[v];
    faceCount = totvNum / 3;
    edgeCount = (totvNum + bdryCount) / 2;

    // alpha resolution -- identical logic to readpack.m's, now shared by
    // both readpack() (which reaches here via the FLOWERS:/BOUQUET: branch
    // above) and loadComplex().
    if (alpha == 0 || (alpha >= 1 && alpha <= nodeCount && bdryFlags[alpha] != 0)) {
        Index candidate = 0;
        for (Index v = 1; v <= nodeCount; ++v) {
            if (bdryFlags[v] == 0) { candidate = -v; break; }
        }
        alpha = candidate;
    }

    if (alpha == 0) {
        std::fprintf(stderr, "Error, no interior vertex was found, stop processing.\n");
        return 0;
    } else if (alpha < 0) {
        if (bdryCount == 0) {
            alpha = 1;
        } else {
            std::vector<Index> seeds;
            for (Index j = 1; j <= nodeCount; ++j) {
                if (bdryFlags[j] == 1) seeds.push_back(j);
            }
            Index a = farVert(seeds);
            alpha = (a < 0) ? 1 : a;
        }
    }

    return nodeCount;
}

Index Packer::finalizeComplex(const std::vector<Scalar>* initRadii,
                               const std::vector<Complex>* initCenters,
                               const std::vector<Scalar>* vAimsIn, const std::string& label) {
    complexCount();
    if (orphanCount > 0) {
        std::fprintf(stderr,
                     "Warning: orphan vertices were found (vertices w/o interior "
                     "neighbors)\n");
    }
    centers.assign(static_cast<size_t>(nodeCount) + 1, Complex(0.0, 0.0));
    radii.assign(static_cast<size_t>(nodeCount) + 1, 0.5);
    if (initRadii) origRadii = *initRadii;
    if (initCenters) origCenters = *initCenters;

    // Mirrors readpack()'s own geometry-dependent handling exactly: given
    // radii/centers are applied directly for Euclidean input, converted via
    // hToEData for Hyperbolic input (needs both radii AND centers to
    // convert), and -- matching the original behavior, not an oversight --
    // simply NOT applied for Spherical input, which falls back to the
    // defaults above even when radii/centers are given.
    if (hes == Geometry::Euclidean) {
        if (initRadii) radii = *initRadii;
        if (initCenters) centers = *initCenters;
    } else if (hes == Geometry::Hyperbolic && initRadii && initCenters) {
        for (Index v = 1; v <= nodeCount; ++v) {
            auto [ez, er] = geom::hToEData((*initCenters)[v], (*initRadii)[v]);
            radii[v] = er;
            centers[v] = ez;
        }
    }

    localcenters = centers;
    localradii = radii;

    if (vAimsIn) {
        vAims = *vAimsIn;
    } else if (vAims.empty()) {
        vAims.assign(static_cast<size_t>(nodeCount) + 1, kTwoPi);
        for (Index k = 1; k <= nodeCount; ++k) {
            if (bdryFlags[k] != 0) vAims[k] = -1.0;
        }
    }

    mode = 1;
    indxMatrices();

    if (!label.empty()) {
        const char* gem = "Euclidean";
        if (hes == Geometry::Hyperbolic) gem = "Hyperbolic";
        else if (hes == Geometry::Spherical) gem = "Spherical";
        std::fprintf(stdout, "Packing %s (%s) is loaded, max pack mode\n", label.c_str(), gem);
    }

    return nodeCount;
}

Index Packer::loadComplex(Index nodeCountIn, const std::vector<std::vector<Index>>& flowersIn,
                           Geometry geometryIn, Index alphaIn, Index gammaIn,
                           const std::vector<Scalar>* initRadii,
                           const std::vector<Complex>* initCenters,
                           const std::vector<Index>* vlistIn,
                           const std::vector<Scalar>* vAimsIn, const std::string& label) {
    if (nodeCountIn <= 0) {
        std::fprintf(stderr, "Packer::loadComplex: nodeCountIn must be positive (got %d)\n",
                     nodeCountIn);
        return 0;
    }
    if (flowersIn.size() != static_cast<size_t>(nodeCountIn) + 1) {
        std::fprintf(stderr,
                     "Packer::loadComplex: flowersIn.size() must be nodeCountIn+1 "
                     "(index 0 unused) -- got %zu, expected %d\n", flowersIn.size(),
                     nodeCountIn + 1);
        return 0;
    }

    cleanse(); // resets alpha/gamma/hes/vAims/vlist/etc. to defaults; must
               // run before we set the fields below, or it would clobber them.
    hes = geometryIn;
    alpha = alphaIn;
    gamma = gammaIn;
    fileName = label.empty() ? "noname" : label;

    if (ingestFlowers(nodeCountIn, flowersIn) == 0) {
        return 0;
    }

    if (vlistIn) vlist = *vlistIn;

    return finalizeComplex(initRadii, initCenters, vAimsIn, label);
}

Index Packer::parseTriangles() {
    throw NotImplementedError(
        "Packer::parseTriangles: bare triangle-list / OFF input reading has not "
        "been ported yet (parse_triangles.m). readpack() supports the *.p "
        "FLOWERS format directly, which is GOPack's preferred format.");
}

std::vector<Complex> Packer::loadTangency(const std::vector<Complex>& centersIn,
                                           const std::vector<Scalar>& radiiIn) const {
    Index tcount = 0;
    for (Index v = 1; v <= nodeCount; ++v) {
        const auto& flower = flowers[v];
        Index n = static_cast<Index>(flower.size());
        if (flower.front() == flower.back()) n -= 1;
        for (Index j = 0; j < n; ++j) {
            if (flower[static_cast<size_t>(j)] > v) tcount++;
        }
    }

    std::vector<Complex> T;
    T.reserve(static_cast<size_t>(tcount));
    for (Index v = 1; v <= nodeCount; ++v) {
        const auto& flower = flowers[v];
        Index n = static_cast<Index>(flower.size());
        if (flower.front() == flower.back()) n -= 1;
        Complex Z = centersIn[v];
        Scalar R = radiiIn[v];
        for (Index j = 0; j < n; ++j) {
            Index w = flower[static_cast<size_t>(j)];
            if (w > v) {
                Complex W = centersIn[w];
                Scalar s = radiiIn[w] + R;
                T.push_back(Z + (R / s) * (W - Z));
            }
        }
    }
    return T;
}

Index Packer::writepack(const std::string& fname, bool euclFlag) {
    std::ofstream out(fname);
    if (!out) {
        std::fprintf(stderr, "Failed open file %s for packing\n", fname.c_str());
        return 0;
    }

    out << "NODECOUNT: " << nodeCount << "\n";
    out << "GEOMETRY: ";
    if (euclFlag) {
        out << "eucl\n";
    } else if (hes == Geometry::Hyperbolic) {
        out << "hyp\n";
    } else if (hes == Geometry::Spherical) {
        out << "sph\n";
    } else {
        out << "eucl\n";
    }
    if (alpha > 0) {
        out << "ALPHA/BETA/GAMMA: " << alpha << " " << 0 << " " << gamma << "\n";
    }
    out << "FLOWERS:\n";
    for (Index k = 1; k <= nodeCount; ++k) {
        out << k << " " << vNum[k] << "  ";
        const auto& flower = flowers[k];
        for (Index j = 0; j <= vNum[k]; ++j) out << flower[static_cast<size_t>(j)] << " ";
        out << "\n";
    }
    out << "RADII: \n";

    std::vector<Scalar> rad2store = radii;
    std::vector<Complex> cent2store = centers;

    if (!euclFlag) {
        if (hes == Geometry::Spherical) {
            std::vector<Complex> T = loadTangency(centers, radii);
            auto [A, B] = geom::affineNormalizer(T);

            std::vector<Complex> Z(static_cast<size_t>(nodeCount) + 1);
            std::vector<Scalar> R(static_cast<size_t>(nodeCount) + 1);
            for (Index v = 1; v <= nodeCount; ++v) {
                auto [sz, sr] = geom::eToSData(A * centers[v] + B, A * radii[v]);
                Z[v] = sz;
                R[v] = sr;
            }
            cent2store = Z;
            rad2store = R;
        } else if (hes == Geometry::Hyperbolic) {
            for (Index v = 1; v <= nodeCount; ++v) {
                auto [hz, hr] = geom::eToHData(centers[v], radii[v]);
                cent2store[v] = hz;
                rad2store[v] = hr;
            }
        }
    }

    for (Index v = 1; v <= nodeCount; ++v) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%16.10f\n", rad2store[v]);
        out << buf;
    }
    out << "CENTERS:\n";
    for (Index v = 1; v <= nodeCount; ++v) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%16.10f %16.10f\n", cent2store[v].real(),
                     cent2store[v].imag());
        out << buf;
    }

    if (!vlist.empty()) {
        out << "VERT_LIST:\n";
        for (Index v : vlist) out << v << "\n";
    }

    bool ndflag = false;
    for (Index v = 1; v <= nodeCount; ++v) {
        const auto& flower = flowers[v];
        bool isBdry = flower.front() != flower.back();
        if (isBdry && vAims[v] >= 0) {
            if (!ndflag) { out << "ANGLE_AIMS:\n"; ndflag = true; }
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%d %13.6f\n", v, vAims[v]);
            out << buf;
        } else if (!isBdry && std::abs(vAims[v] - kTwoPi) > 0.00001) {
            if (!ndflag) { out << "ANGLE_AIMS:\n"; ndflag = true; }
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%d %f\n", v, vAims[v]);
            out << buf;
        }
    }

    out << "\n\nEND\n";
    out.close();

    const char* gem = "eucl";
    if (hes == Geometry::Hyperbolic && !euclFlag) gem = "hyp";
    else if (hes == Geometry::Spherical && !euclFlag) gem = "sph";
    std::fprintf(stdout, "%s packing written to %s\n", gem, fname.c_str());

    return nodeCount;
}

Index Packer::writeEucl(const std::string& fname) {
    return writepack(fname, true);
}

} // namespace gopack
