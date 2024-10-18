#include "pch.h"

// See this SO post: https://stackoverflow.com/questions/73494386/lnk2001-linker-error-while-linking-google-benchmark-lib
#pragma comment ( lib, "Shlwapi.lib" )
#include <benchmark/benchmark.h>

using namespace myhf;

// Overlap Matrix --------------------------------------------------------------------------
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
static void Overlap_CH4_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_3G);
}
static void Overlap_CH4_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_6G);
}

static void Overlap_HH_STO3G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix_2(atoms, STO_3G);
}
static void Overlap_HH_STO6G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix_2(atoms, STO_6G);
}
static void Overlap_H2O_STO3G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix_2(atoms, STO_3G);
}
static void Overlap_H2O_STO6G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix_2(atoms, STO_6G);
}
static void Overlap_CH4_STO3G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix_2(atoms, STO_3G);
}
static void Overlap_CH4_STO6G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix_2(atoms, STO_6G);
}

// Data driven test
static void Overlap_HH_STO3G_TEST(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = test::OverlapMatrix(atoms);
}
static void Overlap_H2O_STO3G_TEST(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = test::OverlapMatrix(atoms);
}
static void Overlap_H2O2_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_3G);
}
static void Overlap_H2O2_STO3G_TEST(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = test::OverlapMatrix(atoms);
}
static void Overlap_H2O2N_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } },
		{ ATOM_TYPE::Nitrogen, 7,   { -1.2, 0.123, -0.6 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = OverlapMatrix(atoms, STO_3G);
}
static void Overlap_H2O2N_STO3G_TEST(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } },
		{ ATOM_TYPE::Nitrogen, 7,   { -1.2, 0.123, -0.6 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = test::OverlapMatrix(atoms);
}

// Kinetic Matrix --------------------------------------------------------------------------
static void Kinetic_HH_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = KineticEnergyMatrix(atoms, STO_3G);
}
static void Kinetic_HH_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = KineticEnergyMatrix(atoms, STO_6G);
}
static void Kinetic_H2O_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = KineticEnergyMatrix(atoms, STO_3G);
}
static void Kinetic_H2O_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = KineticEnergyMatrix(atoms, STO_6G);
}
static void Kinetic_CH4_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = KineticEnergyMatrix(atoms, STO_3G);
}
static void Kinetic_CH4_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd overlap = KineticEnergyMatrix(atoms, STO_6G);
}

// Overlap AND Kinetic Matrix --------------------------------------------------------------------------
static void Overlap_and_Kinetic_HH_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	Eigen::MatrixXd overlap, kinetic;
	for (auto _ : state)
		OverlapAndKineticEnergyMatrix(atoms, STO_3G, overlap, kinetic);
}
static void Overlap_and_Kinetic_HH_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	Eigen::MatrixXd overlap, kinetic;
	for (auto _ : state)
		OverlapAndKineticEnergyMatrix(atoms, STO_6G, overlap, kinetic);
}
static void Overlap_and_Kinetic_H2O_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	Eigen::MatrixXd overlap, kinetic;
	for (auto _ : state)
		OverlapAndKineticEnergyMatrix(atoms, STO_3G, overlap, kinetic);
}
static void Overlap_and_Kinetic_H2O_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	Eigen::MatrixXd overlap, kinetic;
	for (auto _ : state)
		OverlapAndKineticEnergyMatrix(atoms, STO_6G, overlap, kinetic);
}
static void Overlap_and_Kinetic_CH4_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	Eigen::MatrixXd overlap, kinetic;
	for (auto _ : state)
		OverlapAndKineticEnergyMatrix(atoms, STO_3G, overlap, kinetic);
}
static void Overlap_and_Kinetic_CH4_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	Eigen::MatrixXd overlap, kinetic;
	for (auto _ : state)
		OverlapAndKineticEnergyMatrix(atoms, STO_6G, overlap, kinetic);
}

// Nuclear Electron Attraction Matrix --------------------------------------------------------------------------
static void NuclearElectronAttraction_HH_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_3G);
}
static void NuclearElectronAttraction_HH_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_6G);
}
static void NuclearElectronAttraction_H2O_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_3G);
}
static void NuclearElectronAttraction_H2O_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_6G);
}
static void NuclearElectronAttraction_CH4_STO3G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_3G);
}
static void NuclearElectronAttraction_CH4_STO6G(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_6G);
}

static void NuclearElectronAttraction_HH_STO3G_Par(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_Par(atoms, STO_3G);
}
static void NuclearElectronAttraction_HH_STO6G_Par(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_Par(atoms, STO_6G);
}
static void NuclearElectronAttraction_H2O_STO3G_Par(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_Par(atoms, STO_3G);
}
static void NuclearElectronAttraction_H2O_STO6G_Par(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_Par(atoms, STO_6G);
}
static void NuclearElectronAttraction_CH4_STO3G_Par(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_Par(atoms, STO_3G);
}
static void NuclearElectronAttraction_CH4_STO6G_Par(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_Par(atoms, STO_6G);
}





static void NuclearElectronAttraction_HH_STO3G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_3G);
}
static void NuclearElectronAttraction_HH_STO6G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_6G);
}
static void NuclearElectronAttraction_H2O_STO3G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_3G);
}
static void NuclearElectronAttraction_H2O_STO6G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_6G);
}
static void NuclearElectronAttraction_CH4_STO3G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_3G);
}
static void NuclearElectronAttraction_CH4_STO6G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_6G);
}
static void NuclearElectronAttraction_H2O2N_STO3G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } },
		{ ATOM_TYPE::Nitrogen, 7,   { -1.2, 0.123, -0.6 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_3G);
}
static void NuclearElectronAttraction_H2O2N_STO6G_2(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } },
		{ ATOM_TYPE::Nitrogen, 7,   { -1.2, 0.123, -0.6 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_6G);
}

static void NuclearElectronAttraction_HH_STO3G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_3G);
}
static void NuclearElectronAttraction_HH_STO6G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
		{ ATOM_TYPE::Hydrogen, 1, { 1.0, 0.0, 0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_6G);
}
static void NuclearElectronAttraction_H2O_STO3G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_3G);
}
static void NuclearElectronAttraction_H2O_STO6G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, 0.0, 0.24026010 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_6G);
}
static void NuclearElectronAttraction_CH4_STO3G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_3G);
}
static void NuclearElectronAttraction_CH4_STO6G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276,  0.6276,  0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, {  0.6276, -0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276,  0.6276, -0.6276 } },
		{ ATOM_TYPE::Hydrogen, 1, { -0.6276, -0.6276,  0.6276 } },
		{ ATOM_TYPE::Carbon, 6,   {  0.0,     0.0,     0.0 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_6G);
}
static void NuclearElectronAttraction_H2O2N_STO3G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } },
		{ ATOM_TYPE::Nitrogen, 7,   { -1.2, 0.123, -0.6 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_3G);
}
static void NuclearElectronAttraction_H2O2N_STO6G_3(benchmark::State& state)
{
	std::vector<Atom> atoms
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8,   { 0.0, -0.5, 0.24026010 } },
		{ ATOM_TYPE::Oxygen, 8,   { 1.2, 0.0, 0.24026010 } },
		{ ATOM_TYPE::Nitrogen, 7,   { -1.2, 0.123, -0.6 } }
	};

	for (auto _ : state)
		Eigen::MatrixXd neMatrix = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_6G);
}





BENCHMARK(NuclearElectronAttraction_HH_STO3G_2);
BENCHMARK(NuclearElectronAttraction_HH_STO3G_3);

BENCHMARK(NuclearElectronAttraction_HH_STO6G_2);
BENCHMARK(NuclearElectronAttraction_HH_STO6G_3);

BENCHMARK(NuclearElectronAttraction_H2O_STO3G_2);
BENCHMARK(NuclearElectronAttraction_H2O_STO3G_3);

BENCHMARK(NuclearElectronAttraction_H2O_STO6G_2);
BENCHMARK(NuclearElectronAttraction_H2O_STO6G_3);

BENCHMARK(NuclearElectronAttraction_CH4_STO3G_2);
BENCHMARK(NuclearElectronAttraction_CH4_STO3G_3);

BENCHMARK(NuclearElectronAttraction_CH4_STO6G_2);
BENCHMARK(NuclearElectronAttraction_CH4_STO6G_3);

BENCHMARK(NuclearElectronAttraction_H2O2N_STO3G_2);
BENCHMARK(NuclearElectronAttraction_H2O2N_STO3G_3);

BENCHMARK(NuclearElectronAttraction_H2O2N_STO6G_2);
BENCHMARK(NuclearElectronAttraction_H2O2N_STO6G_3);


//// Nuclear Electron Attraction - STO_3G
//BENCHMARK(NuclearElectronAttraction_HH_STO3G);
//BENCHMARK(NuclearElectronAttraction_HH_STO3G_Par);
//
//BENCHMARK(NuclearElectronAttraction_H2O_STO3G);
//BENCHMARK(NuclearElectronAttraction_H2O_STO3G_Par);
//
//BENCHMARK(NuclearElectronAttraction_CH4_STO3G);
//BENCHMARK(NuclearElectronAttraction_CH4_STO3G_Par);
//
//// Nuclear Electron Attraction - STO_6G
//BENCHMARK(NuclearElectronAttraction_HH_STO6G);
//BENCHMARK(NuclearElectronAttraction_HH_STO6G_Par);
//
//BENCHMARK(NuclearElectronAttraction_H2O_STO6G);
//BENCHMARK(NuclearElectronAttraction_H2O_STO6G_Par);
//
//BENCHMARK(NuclearElectronAttraction_CH4_STO6G);
//BENCHMARK(NuclearElectronAttraction_CH4_STO6G_Par);



//// Overlap - STO_3G
//BENCHMARK(Overlap_HH_STO3G);
//BENCHMARK(Overlap_H2O_STO3G);
//BENCHMARK(Overlap_CH4_STO3G);
//
//// Overlap - STO_6G
//BENCHMARK(Overlap_HH_STO6G);
//BENCHMARK(Overlap_H2O_STO6G);
//BENCHMARK(Overlap_CH4_STO6G);
//
//// Kinetic - STO_3G
//BENCHMARK(Kinetic_HH_STO3G);
//BENCHMARK(Kinetic_H2O_STO3G);
//BENCHMARK(Kinetic_CH4_STO3G);
//
//// Kinetic - STO_6G
//BENCHMARK(Kinetic_HH_STO6G);
//BENCHMARK(Kinetic_H2O_STO6G);
//BENCHMARK(Kinetic_CH4_STO6G);
//
//// Overlap & Kinetic - STO_3G
//BENCHMARK(Overlap_and_Kinetic_HH_STO3G);
//BENCHMARK(Overlap_and_Kinetic_H2O_STO3G);
//BENCHMARK(Overlap_and_Kinetic_CH4_STO3G);
//
//// Overlap & Kinetic - STO_6G
//BENCHMARK(Overlap_and_Kinetic_HH_STO6G);
//BENCHMARK(Overlap_and_Kinetic_H2O_STO6G);
//BENCHMARK(Overlap_and_Kinetic_CH4_STO6G);
//
//// Nuclear Electron Attraction - STO_3G
//BENCHMARK(NuclearElectronAttraction_HH_STO3G);
//BENCHMARK(NuclearElectronAttraction_H2O_STO3G);
//BENCHMARK(NuclearElectronAttraction_CH4_STO3G);
//
//// Nuclear Electron Attraction - STO_6G
//BENCHMARK(NuclearElectronAttraction_HH_STO6G);
//BENCHMARK(NuclearElectronAttraction_H2O_STO6G);
//BENCHMARK(NuclearElectronAttraction_CH4_STO6G);


BENCHMARK_MAIN();