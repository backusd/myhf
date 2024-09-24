#include "pch.h"
#include "MatrixTests.h"
#include "Log.h"

using namespace myhf;
using json = nlohmann::json;

using Eigen::MatrixXd;

#define MATRIX_VALUE_TEST(expression, msg) if (!(expression)) { MatrixTestFailed(msg, testNumber, numberOfTests, numSignificantDigits, basis, atoms, expectedMatrix, actualMatrix); ++numFailures; failure = true; }




struct MatrixExpectedDataResult
{
	std::vector<std::string> atoms;
	std::vector<double> positions;
	std::string basis;
	unsigned int numRows;
	unsigned int numCols;
	std::vector<double> values;
};
struct MatrixExpectedData
{
	std::vector<MatrixExpectedDataResult> results;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MatrixExpectedDataResult, atoms, positions, basis, numRows, numCols, values)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(MatrixExpectedData, results)



static bool EqualUpToNSignificantDigits(double v1, double v2, unsigned int numDigits = 6)
{
	// See: https://stackoverflow.com/questions/17380970/how-can-i-check-whether-two-numbers-are-within-x-significant-figures-of-the-pr

	// If both values are significantly close to 0, then just return true
	if (v1 < 0.00000000001 && v2 < 0.00000000001)
		return true;

	bool val = std::abs(v2 - v1) < std::pow(0.1, numDigits) * std::max(std::abs(v1), std::abs(v2));

	if (!val)
	{
		double n1 = std::abs(v2 - v1);
		double n2 = std::pow(0.1, numDigits);
		double n3 = std::max(std::abs(v1), std::abs(v2));
		double n4 = n2 * n3;
		int ijij = 0;
	}

	return val;
}
static bool WithinMargin(double v1, double v2, double margin = 0.0001)
{
	if (v1 == v2) 
		return true;

	if (std::abs(v1) < 1e-10 && std::abs(v2) < 1e-10)
		return true;

	return std::abs(v2 - v1) / std::max(std::abs(v1), std::abs(v2)) < margin;
}
//static constexpr bool WithinMargin(double actual, double expected, double epsilon = 0.0000001) noexcept
//{
//	return actual + epsilon >= expected && actual - epsilon <= expected;
//}

static void MatrixTestFailed(std::string_view msg, unsigned int testNumber, unsigned int numberOfTests, unsigned int numSignificantDigits, const Basis& basis, std::span<Atom> atoms, const MatrixXd& expected, const MatrixXd& actual) noexcept
{
	LOG_ERROR("[FAILED] Test {0} / {1}: {2}", testNumber, numberOfTests, msg);
	LOG_ERROR("Basis: {0}", basis.name);
	LOG_ERROR("{0}", "Atoms:");
	for (const Atom& atom : atoms)
		LOG_ERROR("\t{0} - {1}", ToStringShortName(atom.type), atom.position);
	LOG_ERROR("Required Number of Significant Digits: {0}", numSignificantDigits);
	LOG_ERROR("Expected Matrix:\n\n{0}\n", expected);
	LOG_ERROR("Computed Matrix:\n\n{0}", actual);
}

static int OverlapMatrixTest(const std::string& file)
{
	MatrixExpectedData expectedData;
	try
	{
		std::println("Loading test descriptions: {0}...", file);
		std::ifstream f(file);
		json data = json::parse(f);
		expectedData = data.template get<MatrixExpectedData>();
	}
	catch (const std::exception& e)
	{
		std::println("Failed to read {0}. Caught exception with message: {1}", file, e.what());
		return 1;
	}

	Eigen::MatrixXd expectedMatrix;
//	const double epsilon = 0.000005;
	unsigned int numSignificantDigits = 9;

	const unsigned int numberOfTests = static_cast<unsigned int>(expectedData.results.size());
	unsigned int testNumber = 0;
	unsigned int numFailures = 0;

	for (auto& result : expectedData.results)
	{
		if (numFailures > 0)
			break;

		++testNumber;

		// Construct the overlap matrix
		expectedMatrix = Eigen::MatrixXd(result.numRows, result.numCols);
		for (unsigned int row = 0; row < result.numRows; ++row)
			for (unsigned int col = 0; col < result.numCols; ++col)
				expectedMatrix(row, col) = result.values[static_cast<size_t>(row) * result.numCols + col];

		try
		{
			bool failure = false;

			std::vector<Atom> atoms;
			unsigned int atomCount = static_cast<unsigned int>(result.atoms.size());
			atoms.reserve(atomCount);

			for (unsigned int iii = 0; iii < atomCount; ++iii)
			{
				ATOM_TYPE type = GetAtomType(result.atoms[iii]);
				unsigned int positionIndex0 = iii * 3;
				atoms.emplace_back(type, static_cast<unsigned int>(type), Vec3d{ result.positions[positionIndex0], result.positions[positionIndex0 + 1], result.positions[positionIndex0 + 2] });
			}

			const Basis& basis = GetBasis(result.basis);
			Eigen::MatrixXd actualMatrix = OverlapMatrix(atoms, basis);

			// 1. matrix shape must match
			MATRIX_VALUE_TEST(static_cast<unsigned int>(actualMatrix.cols()) == result.numCols, "Column counts do not match");
			if (failure) continue;

			MATRIX_VALUE_TEST(static_cast<unsigned int>(actualMatrix.rows()) == result.numRows, "Row counts do not match");
			if (failure) continue;

			// 2. Each value must be within the margin for error
			for (unsigned int row_i = 0; row_i < result.numRows; ++row_i)
			{
				for (unsigned int col_i = 0; col_i < result.numCols; ++col_i)
				{
					MATRIX_VALUE_TEST(EqualUpToNSignificantDigits(actualMatrix(row_i, col_i), expectedMatrix(row_i, col_i), numSignificantDigits), std::format("Overlap values at ({0}, {1}) do not match", row_i, col_i));
					if (failure)
						break;
				}
				if (failure) break;
			}
			if (failure) continue;

			// Report Success
			std::string allNames;
			for (const std::string& shortName : result.atoms)
			{
				allNames += shortName;
				allNames += ' ';
			}
			LOG_INFO("[PASSED] Test {0} / {1}: [{2}] {3}", testNumber, numberOfTests, basis.name, allNames);
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

int OverlapMatrixTest()
{
	std::println("==============================================================================");
	std::println("Overlap Matrix Test");
	std::println("==============================================================================");

	int failures = OverlapMatrixTest("data/overlap-matrix-expected-values_two-atom_sto-3g.json");
	failures += OverlapMatrixTest("data/overlap-matrix-expected-values_two-atom_sto-6g.json");
	failures += OverlapMatrixTest("data/overlap-matrix-expected-values_three-atom_sto-3g.json");
	failures += OverlapMatrixTest("data/overlap-matrix-expected-values_three-atom_sto-6g.json");

	return failures;
}

static int KineticEnergyMatrixTest(const std::string& file)
{
	MatrixExpectedData expectedData;
	try
	{
		std::println("Loading test descriptions: {0}...", file);
		std::ifstream f(file);
		json data = json::parse(f);
		expectedData = data.template get<MatrixExpectedData>();
	}
	catch (const std::exception& e)
	{
		std::println("Failed to read {0}. Caught exception with message: {1}", file, e.what());
		return 1;
	}

	Eigen::MatrixXd expectedMatrix;
//	const double epsilon = 0.000005;
	unsigned int numSignificantDigits = 9;

	const unsigned int numberOfTests = static_cast<unsigned int>(expectedData.results.size());
	unsigned int testNumber = 0;
	unsigned int numFailures = 0;

	for (auto& result : expectedData.results)
	{
		if (numFailures > 0)
			break;

		++testNumber;

		// Construct the overlap matrix
		expectedMatrix = Eigen::MatrixXd(result.numRows, result.numCols);
		for (unsigned int row = 0; row < result.numRows; ++row)
			for (unsigned int col = 0; col < result.numCols; ++col)
				expectedMatrix(row, col) = result.values[static_cast<size_t>(row) * result.numCols + col];

		try
		{
			bool failure = false;

			std::vector<Atom> atoms;
			unsigned int atomCount = static_cast<unsigned int>(result.atoms.size());
			atoms.reserve(atomCount);

			for (unsigned int iii = 0; iii < atomCount; ++iii)
			{
				ATOM_TYPE type = GetAtomType(result.atoms[iii]);
				unsigned int positionIndex0 = iii * 3;
				atoms.emplace_back(type, static_cast<unsigned int>(type), Vec3d{ result.positions[positionIndex0], result.positions[positionIndex0 + 1], result.positions[positionIndex0 + 2] });
			}

			const Basis& basis = GetBasis(result.basis);
			Eigen::MatrixXd actualMatrix = KineticEnergyMatrix(atoms, basis);

			// 1. matrix shape must match
			MATRIX_VALUE_TEST(static_cast<unsigned int>(actualMatrix.cols()) == result.numCols, "Column counts do not match");
			if (failure) continue;

			MATRIX_VALUE_TEST(static_cast<unsigned int>(actualMatrix.rows()) == result.numRows, "Row counts do not match");
			if (failure) continue;

			// 2. Each value must be within the margin for error
			for (unsigned int row_i = 0; row_i < result.numRows; ++row_i)
			{
				for (unsigned int col_i = 0; col_i < result.numCols; ++col_i)
				{
					MATRIX_VALUE_TEST(EqualUpToNSignificantDigits(actualMatrix(row_i, col_i), expectedMatrix(row_i, col_i), numSignificantDigits), std::format("Kinetic energy values at ({0}, {1}) do not match", row_i, col_i));
					if (failure)
						break;
				}
				if (failure) break;
			}
			if (failure) continue;

			// Report Success
			std::string allNames;
			for (const std::string& shortName : result.atoms)
			{
				allNames += shortName;
				allNames += ' ';
			}
			LOG_INFO("[PASSED] Test {0} / {1}: [{2}] {3}", testNumber, numberOfTests, basis.name, allNames);
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

int KineticEnergyMatrixTest()
{
	std::println("==============================================================================");
	std::println("Kinetic Matrix Test");
	std::println("==============================================================================");

	int failures = KineticEnergyMatrixTest("data/kinetic-matrix-expected-values_two-atom_sto-3g.json");
	failures += KineticEnergyMatrixTest("data/kinetic-matrix-expected-values_two-atom_sto-6g.json");
	failures += KineticEnergyMatrixTest("data/kinetic-matrix-expected-values_three-atom_sto-3g.json");
	failures += KineticEnergyMatrixTest("data/kinetic-matrix-expected-values_three-atom_sto-6g.json");

	return failures;
}