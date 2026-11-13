#include "../diagnostics/MPILaunchConfig.h"
#include <cstdio>
#include <iostream>

static void check(bool condition) { if (!condition) throw std::runtime_error("MPI allocation check failed"); }
static void reject(int actual, const char* expected, const char* nodefile)
{
    bool failed = false;
    try { GoldenModesMPI::validateRanks(actual, expected, nodefile); }
    catch (const std::runtime_error&) { failed = true; }
    check(failed);
}
int main()
{
    const char* filename = "mpi_launch_test_nodefile.txt";
    { std::ofstream out(filename); for (int i=0;i<35;++i) out << "same-node\n"; out << " \t\r\n"; }
    check(GoldenModesMPI::expectedRanks(nullptr, filename)==35);
    GoldenModesMPI::validateRanks(35, nullptr, filename);
    GoldenModesMPI::validateRanks(69, "69", nullptr);
    GoldenModesMPI::validateRanks(1, nullptr, nullptr);
    GoldenModesMPI::validateRanks(3, "3", filename); // Intentional subset override.
    reject(1, nullptr, filename); // The observed one-node/35-slot failure.
    reject(1, "69", nullptr);
    for (const char* value : {"0", "-1", "1x", "1.0", "999999999999999999999999"}) reject(1, value, nullptr);
    reject(1, nullptr, "missing-nodefile");
    { std::ofstream out(filename); }
    reject(1, nullptr, filename);
    std::remove(filename);
    std::cout << "MPI CPU-slot counting, allocation mismatches and overrides passed\n";
}
