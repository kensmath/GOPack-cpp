#include "gopack/Packer.h"

#include <cctype>
#include <cstdio>
#include <fstream>
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

} // namespace

Index Packer::readpack(const std::string& fname) {
    std::ifstream in(fname);
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

        if (startsWith(tok, "ALPH")) { // ALPHA/BETA/GAMMA
            Index a = 0, b = 0, c = 0;
            in >> a >> b >> c;
            alpha = a;
            gamma = c;
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
            flowers.assign(static_cast<size_t>(nodeCount) + 1, {});
            vNum.assign(static_cast<size_t>(nodeCount) + 1, 0);
            bdryFlags.assign(static_cast<size_t>(nodeCount) + 1, 0);
            bdryCount = 0;

            std::string restOfLine;
            std::getline(in, restOfLine); // consume rest of "FLOWERS:" line

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
                Index m = static_cast<Index>(nums[1]); // vNum
                vNum[j] = m;
                std::vector<Index> flower;
                flower.reserve(static_cast<size_t>(m) + 1);
                for (size_t idx = 2; idx < nums.size(); ++idx) {
                    flower.push_back(static_cast<Index>(nums[idx]));
                }
                if (flower.front() != flower.back()) { // bdry vertex
                    bdryFlags[j] = 1;
                    bdryCount++;
                }
                flowers[j] = std::move(flower);
            }

            intCount = nodeCount - bdryCount;
            Index totvNum = 0;
            for (Index v = 1; v <= nodeCount; ++v) totvNum += vNum[v];
            faceCount = totvNum / 3;
            edgeCount = (totvNum + bdryCount) / 2;

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
        complexCount();
        if (orphanCount > 0) {
            std::fprintf(stderr,
                         "Warning: orphan vertices were found (vertices w/o interior "
                         "neighbors)\n");
        }
        centers.assign(static_cast<size_t>(nodeCount) + 1, Complex(0.0, 0.0));
        radii.assign(static_cast<size_t>(nodeCount) + 1, 0.5);
        if (haveNewRadii) origRadii = newRadii;
        if (haveNewCenters) origCenters = newCenters;
    }

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

    if (vAims.empty()) {
        vAims.assign(static_cast<size_t>(nodeCount) + 1, kTwoPi);
        for (Index k = 1; k <= nodeCount; ++k) {
            if (bdryFlags[k] != 0) vAims[k] = -1.0;
        }
    }

    mode = 1;
    indxMatrices();

    const char* gem = "Euclidean";
    if (hes == Geometry::Hyperbolic) gem = "Hyperbolic";
    else if (hes == Geometry::Spherical) gem = "Spherical";
    std::fprintf(stdout, "Packing %s (%s) is loaded, max pack mode\n", fileName.c_str(), gem);

    return nodeCount;
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
