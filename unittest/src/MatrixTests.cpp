#include "pch.h"
#include "MatrixTests.h"

#include <myhf.h>
using namespace myhf;

#include <nlohmann/json.hpp>
using json = nlohmann::json;

int OverlapMatrixTest()
{
	std::println("==============================================================================");
	std::println("Overlap Matrix Test");
	std::println("==============================================================================");

	json results = {};
	try
	{
		std::ifstream f("data/overlap-matrix-expected-values.json");
		json data = json::parse(f);
		if (!data.contains("results"))
		{
			std::println("Failed to run the test. data/overlap-matrix-expected-values.json does not contain a 'results' field");
			return 1;
		}
		results = data["results"];

		if (!results.is_array())
		{
			std::println("Failed to run the test. 'results' field in data/overlap-matrix-expected-values.json is not an array");
			return 1;
		}
	}
	catch (const std::exception& e)
	{
		std::println("Failed to read data/overlap-matrix-expected-values.json. Caught exception with message: {}", e.what());
		return 1;
	}

	//std::cout << results.dump(4) << std::endl; 

	int loopCounter = 0;

	for (auto& result : results) 
	{
		std::vector<std::string> atoms;
		std::vector<Vec3d> positions;
		std::string basis;
		unsigned int rowCount = 0;
		unsigned int colCount = 0;
		Eigen::MatrixXd expectedOverlap;

		// Atoms
		if (!result.contains("atoms"))
		{
			std::println("Incomplete result - json data is missing 'atoms' field");
			return 1;
		}
		json j_atoms = result["atoms"];
		if (!j_atoms.is_array())
		{
			std::println("Invalid result - 'atoms' field in the json data is not an array");
			return 1;
		}

		for (auto& atom : j_atoms)
		{
			if (!atom.is_string())
			{
				std::println("Invalid result - An item in the 'atoms' field is not a string");
				return 1;
			}

			atoms.push_back(atom.get<std::string>());
		}

		// Positions
		std::vector<double> _positions;
		if (!result.contains("positions"))
		{
			std::println("Incomplete result - json data is missing 'positions' field");
			return 1;
		}
		json j_positions = result["positions"];
		if (!j_positions.is_array())
		{
			std::println("Invalid result - 'positions' field in the json data is not an array");
			return 1;
		}

		for (auto& pos : j_positions)
		{
			if (!pos.is_number())
			{
				std::println("Invalid result - An item in the 'positions' field is not a float");
				std::println("    json['positions']: {}", j_positions.dump(4));
				return 1;
			}

			_positions.push_back(pos.get<double>());
		}

		if (atoms.size() * 3 != _positions.size())
		{
			std::println("Invalid result - The number of 'positions' values is {}, but should be {} because it should be 3 times the number of atoms ({})", _positions.size(), atoms.size() * 3, atoms.size());
			return 1;
		}

		for (int iii = 0; iii < atoms.size(); ++iii)
			positions.emplace_back(_positions[iii], _positions[iii + 1], _positions[iii + 2]);
		
		// Basis
		if (!result.contains("basis"))
		{
			std::println("Incomplete result - json data is missing 'basis' field");
			return 1;
		}
		if (!result["basis"].is_string())
		{
			std::println("Invalid result - 'basis' field in the json data is not a string");
			return 1;
		}
		basis = result["basis"].get<std::string>();

		// Rows/Columns
		if (!result.contains("overlap_rows"))
		{
			std::println("Incomplete result - json data is missing 'overlap_rows' field");
			return 1;
		}
		if (!result["overlap_rows"].is_number_unsigned())
		{
			std::println("Invalid result - 'overlap_rows' field in the json data is not an unsigned int");
			return 1;
		}
		rowCount = result["overlap_rows"].get<unsigned int>();

		if (!result.contains("overlap_cols"))
		{
			std::println("Incomplete result - json data is missing 'overlap_cols' field");
			return 1;
		}
		if (!result["overlap_cols"].is_number_unsigned())
		{
			std::println("Invalid result - 'overlap_cols' field in the json data is not an unsigned int");
			return 1;
		}
		colCount = result["overlap_cols"].get<unsigned int>();

		if (rowCount == 0)
		{
			std::println("Invalid result - 'overlap_rows' field in the json data is 0");
			return 1;
		}
		if (colCount == 0)
		{
			std::println("Invalid result - 'overlap_cols' field in the json data is 0");
			return 1;
		}
		if (rowCount != colCount)
		{
			std::println("Invalid result - 'overlap_rows' field value ({}) not equal to 'overlap_cols' field value ({})", rowCount, colCount);
			return 1;
		}

		// Overlap Matrix
		if (!result.contains("overlap"))
		{
			std::println("Incomplete result - json data is missing 'overlap' field");
			return 1;
		}
		if (!result["overlap"].is_array())
		{
			std::println("Invalid result - 'overlap' field in the json data is not an array");
			return 1;
		}

		std::vector<double> overlapValues;
		overlapValues.reserve(rowCount * colCount);
		for (auto& o : result["overlap"])
		{
			if (!o.is_number())
			{
				std::println("Invalid result - An item in the 'overlap' field is not a float");
				std::println("    json['overlap']: {}", result["overlap"].dump(4));
				return 1;
			}

			overlapValues.push_back(o.get<double>());
		}

		if (overlapValues.size() != rowCount * colCount)
		{
			std::println("Invalid result - 'overlap' field had {} values but the row/col values were {}/{} so there should have been {} values.", overlapValues.size(), rowCount, colCount, rowCount * colCount);
			return 1;
		}

		

		expectedOverlap = Eigen::MatrixXd(rowCount, colCount);
		for (unsigned int row = 0; row < rowCount; ++row)
			for (unsigned int col = 0; col < colCount; ++col)
				expectedOverlap(row, col) = overlapValues[row * colCount + col];

		std::cout << '[' << expectedOverlap << "]\n";


		++loopCounter;
		if (loopCounter == 6)
			break;
	}



	return 0;

	//std::cout << "Unit test\n";
	//
	//std::array<Atom, 2> atoms =
	//{
	//	STO_3G.GetAtom(ATOM_TYPE::Hydrogen),
	//	STO_3G.GetAtom(ATOM_TYPE::Hydrogen)
	//};
	//atoms[1].position = { 0.5356598775430879, 0, 0 }; //  0.283459 * 2
	//Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms);
	//
	//
	//std::print("{}", overlapMatrix);
}