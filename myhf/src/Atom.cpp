#include "pch.h"
#include "Atom.h"


namespace myhf
{
ATOM_TYPE GetAtomType(std::string_view type)
{
	static constexpr std::array names = { "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
		"Carbon", "Nitrogen", "Oxygen", "Flourine", "Neon", "Sodium", "Magnesium", "Aluminum",
		"Silicon", "Phosphorous", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium" };

	static constexpr std::array shortNames = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", 
		"Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca" };

	if (type.length() == 0)
		throw std::invalid_argument(std::format("{}:{} - Atom type parameter must be non-zero in length", __FILE__, __LINE__));

	if (type.length() <= 2)
	{
		for (unsigned int iii = 0; iii < shortNames.size(); ++iii)
		{
			if (type == shortNames[iii])
				return static_cast<ATOM_TYPE>(iii + 1); 
		}
	}
	else
	{
		for (unsigned int iii = 0; iii < names.size(); ++iii)
		{
			if (type == names[iii])
				return static_cast<ATOM_TYPE>(iii + 1);
		}
	}

	throw std::invalid_argument(std::format("{}:{} - Atom type parameter is unrecognized: {}", __FILE__, __LINE__, type));
}

const char* ToStringShortName(ATOM_TYPE type)
{
	static constexpr std::array names = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
		"Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca" };
	return names[static_cast<unsigned int>(type) - 1];
}
const char* ToStringLongName(ATOM_TYPE type)
{
	static constexpr std::array names = { "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
		"Carbon", "Nitrogen", "Oxygen", "Flourine", "Neon", "Sodium", "Magnesium", "Aluminum",
		"Silicon", "Phosphorous", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium" };
	return names[static_cast<unsigned int>(type) - 1];
}
const char* ToString(ATOM_TYPE type) { return ToStringLongName(type); }
}