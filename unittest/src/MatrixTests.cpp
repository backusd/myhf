#include "pch.h"
#include "MatrixTests.h"
#include "Log.h"

using namespace myhf;
using json = nlohmann::json;

using Eigen::MatrixXd;

#define OVERLAP_TEST(expression, msg) if (!(expression)) { OverlapTestFailed(msg, testNumber, numberOfTests, basis, atoms, expectedOverlap, actualOverlap); ++numFailures; failure = true; }




struct OverlapExpectedDataResult
{
	std::vector<std::string> atoms;
	std::vector<double> positions;
	std::string basis;
	unsigned int overlap_rows;
	unsigned int overlap_cols;
	std::vector<double> overlap;
};
struct OverlapExpectedData
{
	std::vector<OverlapExpectedDataResult> results;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(OverlapExpectedDataResult, atoms, positions, basis, overlap_rows, overlap_cols, overlap)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(OverlapExpectedData, results)




static constexpr bool WithinMargin(double actual, double expected, double epsilon = 0.0000001) noexcept
{
	return actual + epsilon >= expected && actual - epsilon <= expected;
}

static void OverlapTestFailed(std::string_view msg, unsigned int testNumber, unsigned int numberOfTests, const Basis& basis, std::span<Atom> atoms, const MatrixXd& expected, const MatrixXd& actual) noexcept
{
	LOG_ERROR("[FAILED] Test {0} / {1}: {2}", testNumber, numberOfTests, msg);
	LOG_ERROR("Basis: {0}", basis.name);
	LOG_ERROR("{0}", "Atoms:");
	for (const Atom& atom : atoms)
		LOG_ERROR("\t{0} - {1}", ToStringShortName(atom.type), atom.position);
	LOG_ERROR("Expected Overlap Matrix:\n\n{0}\n", expected);
	LOG_ERROR("Computed Overlap Matrix:\n\n{0}", actual);
}

int OverlapMatrixTest()
{
	std::println("==============================================================================");
	std::println("Overlap Matrix Test");
	std::println("==============================================================================");

	OverlapExpectedData expectedData;
	try
	{
		std::ifstream f("data/overlap-matrix-expected-values.json");
		json data = json::parse(f);
		expectedData = data.template get<OverlapExpectedData>();
	}
	catch (const std::exception& e)
	{
		std::println("Failed to read data/overlap-matrix-expected-values.json. Caught exception with message: {}", e.what());
		return 1;
	}

	Eigen::MatrixXd expectedOverlap;
	const double epsilon = 0.000001;

	const unsigned int numberOfTests = static_cast<unsigned int>(expectedData.results.size());
	unsigned int testNumber = 0;
	unsigned int numFailures = 0;

	for (auto& result : expectedData.results)
	{
		if (numFailures > 0)
			break;

		++testNumber;	

		// Construct the overlap matrix
		expectedOverlap = Eigen::MatrixXd(result.overlap_rows, result.overlap_cols);
		for (unsigned int row = 0; row < result.overlap_rows; ++row)
			for (unsigned int col = 0; col < result.overlap_cols; ++col)
				expectedOverlap(row, col) = result.overlap[static_cast<size_t>(row) * result.overlap_cols + col];

		try
		{
			bool failure = false;

			std::vector<Atom> atoms;
			unsigned int atomCount = result.atoms.size();
			atoms.reserve(atomCount);

			for (int iii = 0; iii < atomCount; ++iii)
			{
				ATOM_TYPE type = GetAtomType(result.atoms[iii]);
				unsigned int positionIndex0 = iii * 3;
				atoms.emplace_back(type, static_cast<unsigned int>(type), Vec3d{ result.positions[positionIndex0], result.positions[positionIndex0 + 1], result.positions[positionIndex0 + 2] });
			}

			const Basis& basis = GetBasis(result.basis);
			Eigen::MatrixXd actualOverlap = OverlapMatrix(atoms, basis);

			// 1. matrix shape must match
			OVERLAP_TEST(static_cast<unsigned int>(actualOverlap.cols()) == result.overlap_cols, "Column counts do not match");
			if (failure) continue;

			OVERLAP_TEST(static_cast<unsigned int>(actualOverlap.rows()) == result.overlap_rows, "Row counts do not match");
			if (failure) continue;

			// 2. Each value must be within the margin for error
			for (unsigned int row_i = 0; row_i < result.overlap_rows; ++row_i)
			{
				for (unsigned int col_i = 0; col_i < result.overlap_cols; ++col_i)
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