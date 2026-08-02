// Exercises the actual file-format path (readpack -> setMode -> riffle ->
// writepack) end to end on the same hex-flower complex as test_hex_flower,
// this time loaded from a real *.p file, to make sure the FLOWERS parser
// wires everything up the same way manual construction does.
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

#include "gopack/Packer.h"

namespace {
void checkTrue(bool cond, const char* what) {
    if (!cond) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        std::exit(1);
    }
}
} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <path-to-small_disc.p>\n", argv[0]);
        return 2;
    }

    gopack::Packer packer;
    gopack::Index nc = packer.readpack(argv[1]);
    checkTrue(nc == 7, "readpack should report nodeCount == 7");
    checkTrue(packer.intCount == 1, "expected 1 interior vertex");
    checkTrue(packer.bdryCount == 6, "expected 6 boundary vertices");

    int md = packer.setMode(1);
    checkTrue(md == 1, "setMode(1) should succeed");

    gopack::RiffleResult result = packer.riffle(100);
    checkTrue(result.cycles >= 0, "riffle should not report an error");

    for (gopack::Index v = 1; v <= packer.nodeCount; ++v) {
        checkTrue(packer.radii[v] > 0.0, "all radii should be positive after riffle");
    }

    std::string outPath = (std::filesystem::temp_directory_path() / "small_disc_out.p").string();
    gopack::Index written = packer.writepack(outPath);
    checkTrue(written == 7, "writepack should report nodeCount == 7");

    // Read it back in a fresh Packer to make sure our own output is
    // self-consistent (round-trips through readpack without error).
    gopack::Packer reloaded;
    gopack::Index nc2 = reloaded.readpack(outPath);
    checkTrue(nc2 == 7, "re-reading our own writepack() output should also report 7 vertices");

    std::fprintf(stdout, "readpack_roundtrip: OK (%d riffle passes, output at %s)\n",
                 result.cycles, outPath.c_str());
    return 0;
}
