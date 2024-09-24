#include "pch.h"

// See this SO post: https://stackoverflow.com/questions/73494386/lnk2001-linker-error-while-linking-google-benchmark-lib
#pragma comment ( lib, "Shlwapi.lib" )
#include <benchmark/benchmark.h>

using namespace myhf;

static void Overlap_HH_STO3G(benchmark::State& state) 
{
	std::vector<Atom> atoms
    {
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

    for (auto _ : state)
        Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_3G);
}
BENCHMARK(Overlap_HH_STO3G);

static void Overlap_HH_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_6G);
}
BENCHMARK(Overlap_HH_STO6G);

static void Overlap_H2O_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_3G);
}
BENCHMARK(Overlap_H2O_STO3G);

static void Overlap_H2O_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_6G);
}
BENCHMARK(Overlap_H2O_STO6G);





BENCHMARK_MAIN();