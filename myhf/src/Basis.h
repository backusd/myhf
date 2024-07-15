#pragma once


#include "Shell.h"
#include "Atom.h"



namespace myhf
{
class Basis
{
public:
	const std::string_view name;
	std::vector<Atom> atoms;

	constexpr Atom GetAtom(ATOM_TYPE type) const noexcept 
	{ 
		assert(static_cast<unsigned int>(type) > 0);
		assert(static_cast<unsigned int>(type) - 1 < atoms.size());
		return atoms[static_cast<unsigned int>(type) - 1];
	}



	void Print() const noexcept;

};

const Basis STO3G{ "STO3G",
	// vector<Atom>
	{
		// Atom [Hydrogen]
		{
			// atom type, number of electrons, position
			ATOM_TYPE::Hydrogen, 
			1, 
			{ 0, 0, 0 },

			// std::vector<ContractedGaussianShell>
			{
				// ContractedGaussianShell
				{
					// std::vector<ContractedGaussianOrbital> (basis functions)
					{
						// ContractedGaussianOrbital
						{
							// std::vector<PrimitiveGaussian>
							{
								// PrimitiveGaussians (alpha, coefficient)
								{ 3.42525091, 0.15432897 },
								{ 0.62391373, 0.53532814 },
								{ 0.16885540, 0.44463454 }
							},
							// QuantumNumbers (angular momentum
							{ 0, 0, 0 }
						}
					}
				}
			}
		},

		// Atom [Helium]
		{
			// atom type, number of electrons, position
			ATOM_TYPE::Helium, 
			2, 
			{ 0, 0, 0 },

			// std::vector<ContractedGaussianShell>
			{
				// ContractedGaussianShell
				{
					// std::vector<ContractedGaussianOrbital> (basis functions)
					{
						// ContractedGaussianOrbital
						{
							// std::vector<PrimitiveGaussian>
							{
								// PrimitiveGaussians (alpha, coefficient)
								{ 6.36242139, 0.15432897 },
								{ 1.15892300, 0.53532814 },
								{ 0.31364979, 0.44463454 }
							},
							// QuantumNumbers (angular momentum)
							{ 0, 0, 0 }
						}
					}
				}
			}
		},

		// Atom [Lithium]
		{
			// atom type, number of electrons, position
			ATOM_TYPE::Lithium, 
			3, 
			{ 0, 0, 0 },

			// std::vector<ContractedGaussianShell>
			{
				// ContractedGaussianShell (s)
				{
					// std::vector<ContractedGaussianOrbital> (basis functions)
					{
						// ContractedGaussianOrbital
						{
							// std::vector<PrimitiveGaussian>
							{
								// PrimitiveGaussians (alpha, coefficient)
								{ 16.1195750, 0.15432897 },
								{ 2.9362007, 0.53532814 },
								{ 0.7946505, 0.44463454 }
							},
							// QuantumNumbers (angular momentum)
							{ 0, 0, 0 }
						}
					}
				},

				// ContractedGaussianShell (sp)
				{
					// std::vector<ContractedGaussianOrbital> (basis functions)
					{
						// ContractedGaussianOrbitals
						{
							// std::vector<PrimitiveGaussian>
							{
								// PrimitiveGaussians (alpha, coefficient)
								{ 0.6362897, -0.09996723 },
								{ 0.1478601, 0.39951283 },
								{ 0.0480887, 0.70011547 }
							},
							// QuantumNumbers (angular momentum)
							{ 0, 0, 0 }
						},
						{
							// std::vector<PrimitiveGaussian>
							{
								// PrimitiveGaussians (alpha, coefficient)
								{ 0.6362897, 0.15591627 },
								{ 0.1478601, 0.60768372 },
								{ 0.0480887, 0.39195739 }
							},
							// QuantumNumbers (angular momentum)
							{ 1, 0, 0 }
						},
						{
							// std::vector<PrimitiveGaussian>
							{
								// PrimitiveGaussians (alpha, coefficient)
								{ 0.6362897, 0.15591627 },
								{ 0.1478601, 0.60768372 },
								{ 0.0480887, 0.39195739 }
							},
							// QuantumNumbers (angular momentum)
							{ 0, 1, 0 }
						},
						{
							// std::vector<PrimitiveGaussian>
							{
								// PrimitiveGaussians (alpha, coefficient)
								{ 0.6362897, 0.15591627 },
								{ 0.1478601, 0.60768372 },
								{ 0.0480887, 0.39195739 }
							},
							// QuantumNumbers (angular momentum)
							{ 0, 0, 1 }
						},
					}
				}
			}
		}
	}
};


}