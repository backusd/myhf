#include "pch.h"
#include "MatrixTests.h"
#include "Log.h"

#include <myhf.h>
using namespace myhf;

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#define OVERLAP_TEST(expression, msg) if (!(expression)) { OverlapTestFailed(msg, testNumber, numberOfTests, molec, expectedOverlap, actualOverlap); ++numFailures; failure = true; }

static constexpr bool WithinMargin(double actual, double expected, double epsilon = 0.0000001) noexcept
{
	return actual + epsilon >= expected && actual - epsilon <= expected;
}

static void OverlapTestFailed(std::string_view msg, unsigned int testNumber, unsigned int numberOfTests, const Molecule& molec, const Eigen::MatrixXd& expected, const Eigen::MatrixXd& actual) noexcept
{
	LOG_ERROR("[FAILED] Test {0} / {1}: {2}", testNumber, numberOfTests, msg);
	LOG_ERROR("Molecule:\n{0}", molec);
	LOG_ERROR("Expected Overlap Matrix:\n\n{0}\n", expected);
	LOG_ERROR("Computed Overlap Matrix:\n\n{0}", actual);
}

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

	std::vector<std::string> atomNames;
	std::vector<Vec3d> positions;
	std::string basis;
	unsigned int rowCount = 0;
	unsigned int colCount = 0;
	Eigen::MatrixXd expectedOverlap;
	const double epsilon = 0.000001;

	const unsigned int numberOfTests = static_cast<unsigned int>(results.size());
	unsigned int testNumber = 0;
	unsigned int numFailures = 0;

	for (auto& result : results) 
	{
		if (numFailures > 0)
			break;

		++testNumber;

		atomNames.clear();
		positions.clear();
		basis.clear();		

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

			atomNames.push_back(atom.get<std::string>());
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

		if (atomNames.size() * 3 != _positions.size())
		{
			std::println("Invalid result - The number of 'positions' values is {}, but should be {} because it should be 3 times the number of atoms ({})", _positions.size(), atomNames.size() * 3, atomNames.size());
			return 1;
		}

		for (int iii = 0; iii < atomNames.size(); ++iii)
			positions.emplace_back(_positions[iii], _positions[static_cast<size_t>(iii) + 1], _positions[static_cast<size_t>(iii) + 2]);
		
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
		overlapValues.reserve(static_cast<size_t>(rowCount) * colCount);
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

		if (overlapValues.size() != static_cast<size_t>(rowCount) * colCount)
		{
			std::println("Invalid result - 'overlap' field had {} values but the row/col values were {}/{} so there should have been {} values.", overlapValues.size(), rowCount, colCount, rowCount * colCount);
			return 1;
		}

		expectedOverlap = Eigen::MatrixXd(rowCount, colCount);
		for (unsigned int row = 0; row < rowCount; ++row)
			for (unsigned int col = 0; col < colCount; ++col)
				expectedOverlap(row, col) = overlapValues[static_cast<size_t>(row) * colCount + col];


		// Expected Overlap matrix has been parsed... Now do the test
		try
		{
			bool failure = false;

			std::vector<Atom> atoms;
			atoms.reserve(atomNames.size());
			for (int iii = 0; iii < atomNames.size(); ++iii)
			{
				ATOM_TYPE type = GetAtomType(atomNames[iii]);
				atoms.emplace_back(type, static_cast<unsigned int>(type), positions[iii]);
			}

			Molecule molec(std::move(atoms), GetBasis(basis));
			Eigen::MatrixXd actualOverlap = molec.OverlapMatrix();

			// 1. matrix shape must match
			OVERLAP_TEST(static_cast<unsigned int>(actualOverlap.cols()) == colCount, "Column counts do not match");
			if (failure) continue;

			OVERLAP_TEST(static_cast<unsigned int>(actualOverlap.rows()) == rowCount, "Row counts do not match");
			if (failure) continue;

			// 2. Each value must be within the margin for error
			for (unsigned int row_i = 0; row_i < rowCount; ++row_i)
			{
				for (unsigned int col_i = 0; col_i < colCount; ++col_i)
				{
					OVERLAP_TEST(WithinMargin(actualOverlap(row_i, col_i), expectedOverlap(row_i, col_i), epsilon), std::format("Overlap values at ({0}, {1}) do not match", row_i, col_i));
					if (failure) 
						break;
				}
				if (failure) break;
			}
			if (failure) continue;

			// Report Success
			std::string allNames;
			for (const std::string& shortName : atomNames)
			{
				allNames += shortName;
				allNames += ' ';
			}
			LOG_INFO("[PASSED] Test {0} / {1}: {2}", testNumber, numberOfTests, allNames); 
		}
		catch (std::exception& ex)
		{
			LOG_ERROR("ERROR: Caught exception: {}", ex.what());
		}
	}

	if (numFailures > 0)
		LOG_WARN("SUMMARY: Passed = {0}   Failed = {1}", numberOfTests - numFailures, numFailures);
	else
		LOG_INFO("SUMMARY: All tests passed");

	return numFailures;
}