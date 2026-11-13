#ifndef GOLDENMODES_MPI_LAUNCH_CONFIG_H
#define GOLDENMODES_MPI_LAUNCH_CONFIG_H

#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

namespace GoldenModesMPI {
// The launcher can explicitly request a subset for a benchmark. Otherwise a PBS
// allocation means one MPI process per nonempty CPU-slot entry, not per hostname.
inline unsigned expectedRanks(const char* explicitRanks, const char* nodefile)
{
    if (explicitRanks && *explicitRanks) {
        const std::string value(explicitRanks);
        if (value.find_first_not_of("0123456789") != std::string::npos)
            throw std::runtime_error("GOLDENMODES_MPI_EXPECTED_RANKS must be a positive integer.");
        unsigned long n;
        try { n = std::stoul(value); }
        catch (const std::exception&) { throw std::runtime_error("Invalid GOLDENMODES_MPI_EXPECTED_RANKS."); }
        if (!n || n > static_cast<unsigned long>(std::numeric_limits<int>::max()))
            throw std::runtime_error("Invalid GOLDENMODES_MPI_EXPECTED_RANKS.");
        return static_cast<unsigned>(n);
    }
    if (!nodefile || !*nodefile) return 0; // Ordinary local runs need no scheduler.
    std::ifstream input(nodefile);
    if (!input) throw std::runtime_error("Cannot read PBS_NODEFILE to verify the MPI allocation.");
    unsigned count = 0;
    std::string line;
    while (std::getline(input, line)) {
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        if (count == static_cast<unsigned>(std::numeric_limits<int>::max()))
            throw std::runtime_error("Too many CPU slots in PBS_NODEFILE.");
        ++count;
    }
    if (input.bad() || !count) throw std::runtime_error("Invalid or empty PBS_NODEFILE.");
    return count;
}

inline void validateRanks(int actual, const char* explicitRanks, const char* nodefile)
{
    const unsigned expected = expectedRanks(explicitRanks, nodefile);
    if (actual < 1 || (expected && static_cast<unsigned>(actual) != expected))
        throw std::runtime_error("MPI launch mismatch: expected " + std::to_string(expected) +
            " MPI processes, but started " + std::to_string(actual) +
            ". Launch with the CPU-slot count (NCPU), not the node count (NNODES).");
}
}
#endif
