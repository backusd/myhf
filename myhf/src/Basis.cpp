#include "pch.h"
#include "Basis.h"


namespace myhf
{
const Basis& GetBasis(std::string_view basis)
{
	if (basis == "STO-3G" || basis == "sto-3g") return STO_3G;
	if (basis == "STO-6G" || basis == "sto-6g") return STO_6G;

	throw std::invalid_argument(std::format("{}:{} - Basis parameter is unrecognized: {}", __FILE__, __LINE__, basis));
}

}