#pragma once
#include "Shell.h"
#include "Atom.h"

namespace myhf
{
struct BasisAtom
{
	std::vector<ContractedGaussianShell> shells;

	constexpr unsigned int NumberOfContractedGaussians() const noexcept
	{
		unsigned int res = 0;
		for (auto const& shell : shells)
			res += shell.CountNumberOfContractedGaussians();
		return res;
	}
	constexpr QuantumNumbers GetMaxQN(double alpha) const noexcept
	{
		QuantumNumbers maxQN(0, 0, 0);

		for (const auto& shell : shells)
			for (const auto& orbital : shell.basisFunctions)
				for (const auto& gaussian : orbital.gaussianOrbitals)
					if (alpha == gaussian.alpha)
					{
						maxQN.l = std::max(maxQN.l, orbital.angularMomentum.l);
						maxQN.m = std::max(maxQN.m, orbital.angularMomentum.m);
						maxQN.n = std::max(maxQN.m, orbital.angularMomentum.n);
					}

		return maxQN;
	}
};

struct Basis
{
	const std::string_view name;
	std::vector<BasisAtom> atoms;

	constexpr const BasisAtom& GetAtom(ATOM_TYPE type) const noexcept 
	{ 
		assert(static_cast<unsigned int>(type) > 0);
		assert(static_cast<unsigned int>(type) - 1 < atoms.size());
		return atoms[static_cast<unsigned int>(type) - 1];
	}
};

const char* ToString(const Basis& basis);

const Basis& GetBasis(std::string_view basis);

#define BRACES(...) { __VA_ARGS__ }

#define STO_3G_ATOM(atom, details) { details }
#define STO_3G_S(a1, c1, a2, c2, a3, c3) {{{{{a1, c1}, {a2, c2}, {a3, c3}}, {0,0,0} }}}
#define STO_3G_SP(a1, c1_s, c1_p, a2, c2_s, c2_p, a3, c3_s, c3_p) {{ \
	{{{a1, c1_s}, {a2, c2_s}, {a3, c3_s}}, {0,0,0}}, \
	{{{a1, c1_p}, {a2, c2_p}, {a3, c3_p}}, {1,0,0}}, \
	{{{a1, c1_p}, {a2, c2_p}, {a3, c3_p}}, {0,1,0}}, \
	{{{a1, c1_p}, {a2, c2_p}, {a3, c3_p}}, {0,0,1}} }}

const Basis STO_3G{ "STO-3G",
	// std::vector<BasisAtom>
	{
		STO_3G_ATOM(ATOM_TYPE::Hydrogen, BRACES(
			STO_3G_S(
				0.3425250914E+01, 0.1543289673E+00,
				0.6239137298E+00, 0.5353281423E+00,
				0.1688554040E+00, 0.4446345422E+00
			) 
		)),
		STO_3G_ATOM(ATOM_TYPE::Helium, BRACES(
			STO_3G_S(
				0.6362421394E+01, 0.1543289673E+00,
				0.1158922999E+01, 0.5353281423E+00,
				0.3136497915E+00, 0.4446345422E+00
			) 
		)),
		STO_3G_ATOM(ATOM_TYPE::Lithium, BRACES(			
			STO_3G_S(
				0.1611957475E+02, 0.1543289673E+00,
				0.2936200663E+01, 0.5353281423E+00,
				0.7946504870E+00, 0.4446345422E+00
			),
			STO_3G_SP(
				0.6362897469E+00, -0.9996722919E-01, 0.1559162750E+00,
				0.1478600533E+00,  0.3995128261E+00, 0.6076837186E+00,
				0.4808867840E-01,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Beryllium, BRACES(
			STO_3G_S(
				0.3016787069E+02, 0.1543289673E+00,
				0.5495115306E+01, 0.5353281423E+00,
				0.1487192653E+01, 0.4446345422E+00
			),
			STO_3G_SP(
				0.1314833110E+01, -0.9996722919E-01, 0.1559162750E+00,
				0.3055389383E+00,  0.3995128261E+00, 0.6076837186E+00,
				0.9937074560E-01,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Boron, BRACES(
			STO_3G_S(
				0.4879111318E+02, 0.1543289673E+00,
				0.8887362172E+01, 0.5353281423E+00,
				0.2405267040E+01, 0.4446345422E+00
			),
			STO_3G_SP(
				0.2236956142E+01, -0.9996722919E-01, 0.1559162750E+00,
				0.5198204999E+00,  0.3995128261E+00, 0.6076837186E+00,
				0.1690617600E+00,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Carbon, BRACES(
			STO_3G_S(
				0.7161683735E+02, 0.1543289673E+00,
				0.1304509632E+02, 0.5353281423E+00,
				0.3530512160E+01, 0.4446345422E+00
			),
			STO_3G_SP(
				0.2941249355E+01, -0.9996722919E-01, 0.1559162750E+00,
				0.6834830964E+00,  0.3995128261E+00, 0.6076837186E+00,
				0.2222899159E+00,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Nitrogen, BRACES(
			STO_3G_S(
				0.9910616896E+02, 0.1543289673E+00,
				0.1805231239E+02, 0.5353281423E+00,
				0.4885660238E+01, 0.4446345422E+00
			),
			STO_3G_SP(
				0.3780455879E+01, -0.9996722919E-01, 0.1559162750E+00,
				0.8784966449E+00,  0.3995128261E+00, 0.6076837186E+00,
				0.2857143744E+00,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Oxygen, BRACES(
			STO_3G_S(
				0.1307093214E+03, 0.1543289673E+00,
				0.2380886605E+02, 0.5353281423E+00,
				0.6443608313E+01, 0.4446345422E+00
			),
			STO_3G_SP(
				0.5033151319E+01, -0.9996722919E-01, 0.1559162750E+00, 
				0.1169596125E+01,  0.3995128261E+00, 0.6076837186E+00, 
				0.3803889600E+00,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Flourine, BRACES(
			STO_3G_S(
				0.1666791340E+03, 0.1543289673E+00,
				0.3036081233E+02, 0.5353281423E+00,
				0.8216820672E+01, 0.4446345422E+00
			),
			STO_3G_SP(
				0.6464803249E+01, -0.9996722919E-01, 0.1559162750E+00, 
				0.1502281245E+01,  0.3995128261E+00, 0.6076837186E+00, 
				0.4885884864E+00,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Neon, BRACES(
			STO_3G_S(
				0.2070156070E+03, 0.1543289673E+00,
				0.3770815124E+02, 0.5353281423E+00,
				0.1020529731E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.8246315120E+01, -0.9996722919E-01, 0.1559162750E+00,
				0.1916266291E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.6232292721E+00,  0.7001154689E+00, 0.3919573931E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Sodium, BRACES(
			STO_3G_S(
				0.2507724300E+03, 0.1543289673E+00,
				0.4567851117E+02, 0.5353281423E+00,
				0.1236238776E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.1204019274E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.2797881859E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.9099580170E+00,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.1478740622E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.4125648801E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.1614750979E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Magnesium, BRACES(
			STO_3G_S(
				0.2992374137E+03, 0.1543289673E+00,
				0.5450646845E+02, 0.5353281423E+00,
				0.1475157752E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.1512182352E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.3513986579E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.1142857498E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.1395448293E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.3893265318E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.1523797659E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Aluminum, BRACES(
			STO_3G_S(
				0.3514214767E+03, 0.1543289673E+00,
				0.6401186067E+02, 0.5353281423E+00,
				0.1732410761E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.1889939621E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.4391813233E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.1428353970E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.1395448293E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.3893265318E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.1523797659E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Silicon, BRACES(
			STO_3G_S(
				0.4077975514E+03, 0.1543289673E+00,
				0.7428083305E+02, 0.5353281423E+00,
				0.2010329229E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.2319365606E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.5389706871E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.1752899952E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.1478740622E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.4125648801E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.1614750979E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Phosphorous, BRACES(
			STO_3G_S(
				0.4683656378E+03, 0.1543289673E+00,
				0.8531338559E+02, 0.5353281423E+00,
				0.2308913156E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.2803263958E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.6514182577E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.2118614352E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.1743103231E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.4863213771E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.1903428909E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Sulfur, BRACES(
			STO_3G_S(
				0.5331257359E+03, 0.1543289673E+00,
				0.9710951830E+02, 0.5353281423E+00,
				0.2628162542E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.3332975173E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.7745117521E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.2518952599E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.2029194274E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.5661400518E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.2215833792E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Chlorine, BRACES(
			STO_3G_S(
				0.6013456136E+03, 0.1543289673E+00,
				0.1095358542E+03, 0.5353281423E+00,
				0.2964467686E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.3896041889E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.9053563477E+01,  0.3995128261E+00, 0.6076837186E+00,
				0.2944499834E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.2129386495E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.5940934274E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.2325241410E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Argon, BRACES(
			STO_3G_S(
				0.6744465184E+03, 0.1543289673E+00,
				0.1228512753E+03, 0.5353281423E+00,
				0.3324834945E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.4516424392E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.1049519900E+02,  0.3995128261E+00, 0.6076837186E+00,
				0.3413364448E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.2621366518E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.7313546050E+00,  0.2255954336E+00, 0.5951670053E+00,
				0.2862472356E+00,  0.9003984260E+00, 0.4620010120E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Potassium, BRACES(
			STO_3G_S(
				0.7715103681E+03, 0.1543289673E+00,
				0.1405315766E+03, 0.5353281423E+00,
				0.3803332899E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.5240203979E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.1217710710E+02,  0.3995128261E+00, 0.6076837186E+00,
				0.3960373165E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.3651583985E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.1018782663E+01,  0.2255954336E+00, 0.5951670053E+00,
				0.3987446295E+00,  0.9003984260E+00, 0.4620010120E+00
			),
			STO_3G_SP(
				0.5039822505E+00, -0.3088441214E+00, -0.1215468600E+00,
				0.1860011465E+00,  0.1960641165E-01,  0.5715227604E+00,
				0.8214006743E-01,  0.1131034442E+01,  0.5498949471E+00
			)
		)),
		STO_3G_ATOM(ATOM_TYPE::Calcium, BRACES(
			STO_3G_S(
				0.8540324951E+03, 0.1543289673E+00,
				0.1555630851E+03, 0.5353281423E+00,
				0.4210144179E+02, 0.4446345422E+00
			),
			STO_3G_SP(
				0.5956029944E+02, -0.9996722919E-01, 0.1559162750E+00,
				0.1384053270E+02,  0.3995128261E+00, 0.6076837186E+00,
				0.4501370797E+01,  0.7001154689E+00, 0.3919573931E+00
			),
			STO_3G_SP(
				0.4374706256E+01, -0.2196203690E+00, 0.1058760429E-01,
				0.1220531941E+01,  0.2255954336E+00, 0.5951670053E+00,
				0.4777079296E+00,  0.9003984260E+00, 0.4620010120E+00
			),
			STO_3G_SP(
				0.4558489757E+00, -0.3088441214E+00, -0.1215468600E+00,
				0.1682369410E+00,  0.1960641165E-01,  0.5715227604E+00,
				0.7429520696E-01,  0.1131034442E+01,  0.5498949471E+00
			)
		))
	}
};

#define STO_6G_ATOM(atom, details) { details }
#define STO_6G_S(a1, c1, a2, c2, a3, c3, a4, c4, a5, c5, a6, c6) {{{{{a1, c1}, {a2, c2}, {a3, c3}, {a4, c4}, {a5, c5}, {a6, c6}}, {0,0,0} }}}
#define STO_6G_SP(a1, c1_s, c1_p, a2, c2_s, c2_p, a3, c3_s, c3_p, a4, c4_s, c4_p, a5, c5_s, c5_p, a6, c6_s, c6_p) {{ \
	{{{a1, c1_s}, {a2, c2_s}, {a3, c3_s}, {a4, c4_s}, {a5, c5_s}, {a6, c6_s} }, {0,0,0}}, \
	{{{a1, c1_p}, {a2, c2_p}, {a3, c3_p}, {a4, c4_p}, {a5, c5_p}, {a6, c6_p} }, {1,0,0}}, \
	{{{a1, c1_p}, {a2, c2_p}, {a3, c3_p}, {a4, c4_p}, {a5, c5_p}, {a6, c6_p} }, {0,1,0}}, \
	{{{a1, c1_p}, {a2, c2_p}, {a3, c3_p}, {a4, c4_p}, {a5, c5_p}, {a6, c6_p} }, {0,0,1}} }}

const Basis STO_6G{ "STO-6G",
// std::vector<BasisAtom>
{
	STO_6G_ATOM(ATOM_TYPE::Hydrogen, BRACES(
		STO_6G_S(
			0.3552322122E+02, 0.9163596281E-02,
			0.6513143725E+01, 0.4936149294E-01,
			0.1822142904E+01, 0.1685383049E+00,
			0.6259552659E+00, 0.3705627997E+00,
			0.2430767471E+00, 0.4164915298E+00,
			0.1001124280E+00, 0.1303340841E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Helium, BRACES(
		STO_6G_S(
			0.6598456824E+02, 0.9163596281E-02,
			0.1209819836E+02, 0.4936149294E-01,
			0.3384639924E+01, 0.1685383049E+00,
			0.1162715163E+01, 0.3705627997E+00,
			0.4515163224E+00, 0.4164915298E+00,
			0.1859593559E+00, 0.1303340841E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Lithium, BRACES(			
		STO_6G_S(
			0.1671758462E+03, 0.9163596281E-02,
			0.3065150840E+02, 0.4936149294E-01,
			0.8575187477E+01, 0.1685383049E+00,
			0.2945808337E+01, 0.3705627997E+00,
			0.1143943581E+01, 0.4164915298E+00,
			0.4711391391E+00, 0.1303340841E+00
		),
		STO_6G_SP(
			0.6597563981E+01, -0.1325278809E-01, 0.3759696623E-02,
			0.1305830092E+01, -0.4699171014E-01, 0.3767936984E-01,
			0.4058510193E+00, -0.3378537151E-01, 0.1738967435E+00,
			0.1561455158E+00,  0.2502417861E+00, 0.4180364347E+00,
			0.6781410394E-01,  0.5951172526E+00, 0.4258595477E+00,
			0.3108416550E-01,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Beryllium, BRACES(			
		STO_6G_S(
			0.3128704937E+03, 0.9163596281E-02,
			0.5736446253E+02, 0.4936149294E-01,
			0.1604850940E+02, 0.1685383049E+00,
			0.5513096119E+01, 0.3705627997E+00,
			0.2140896553E+01, 0.4164915298E+00,
			0.8817394283E+00, 0.1303340841E+00
		),
		STO_6G_SP(
			0.1363324744E+02, -0.1325278809E-01, 0.3759696623E-02,
			0.2698375464E+01, -0.4699171014E-01, 0.3767936984E-01,
			0.8386530829E+00, -0.3378537151E-01, 0.1738967435E+00,
			0.3226600698E+00,  0.2502417861E+00, 0.4180364347E+00,
			0.1401314882E+00,  0.5951172526E+00, 0.4258595477E+00,
			0.6423251387E-01,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Boron, BRACES(
		STO_6G_S(
			0.5060118369E+03, 0.9163596281E-02,
			0.9277671639E+02, 0.4936149294E-01,
			0.2595558190E+02, 0.1685383049E+00,
			0.8916442908E+01, 0.3705627997E+00,
			0.3462515703E+01, 0.4164915298E+00,
			0.1426055179E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.2319456087E+02, -0.1325278809E-01, 0.3759696623E-02,
			0.4590808918E+01, -0.4699171014E-01, 0.3767936984E-01,
			0.1426819990E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.5489490791E+00,  0.2502417861E+00, 0.4180364347E+00,
			0.2384089591E+00,  0.5951172526E+00, 0.4258595477E+00,
			0.1092802694E+00,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Carbon, BRACES(
		STO_6G_S(
			0.7427370491E+03, 0.9163596281E-02,
			0.1361800249E+03, 0.4936149294E-01,
			0.3809826352E+02, 0.1685383049E+00,
			0.1308778177E+02, 0.3705627997E+00,
			0.5082368648E+01, 0.4164915298E+00,
			0.2093200076E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.3049723950E+02, -0.1325278809E-01, 0.3759696623E-02,
			0.6036199601E+01, -0.4699171014E-01, 0.3767936984E-01,
			0.1876046337E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.7217826470E+00,  0.2502417861E+00, 0.4180364347E+00,
			0.3134706954E+00,  0.5951172526E+00, 0.4258595477E+00,
			0.1436865550E+00,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Nitrogen, BRACES(
		STO_6G_S(
			0.1027828458E+04, 0.9163596281E-02,
			0.1884512226E+03, 0.4936149294E-01,
			0.5272186097E+02, 0.1685383049E+00,
			0.1811138217E+02, 0.3705627997E+00,
			0.7033179691E+01, 0.4164915298E+00,
			0.2896651794E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.3919880787E+02, -0.1325278809E-01, 0.3759696623E-02,
			0.7758467071E+01, -0.4699171014E-01, 0.3767936984E-01,
			0.2411325783E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.9277239437E+00,  0.2502417861E+00, 0.4180364347E+00,
			0.4029111410E+00,  0.5951172526E+00, 0.4258595477E+00,
			0.1846836552E+00,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Oxygen, BRACES(
		STO_6G_S(
			0.1355584234E+04, 0.9163596281E-02,
			0.2485448855E+03, 0.4936149294E-01,
			0.6953390229E+02, 0.1685383049E+00,
			0.2388677211E+02, 0.3705627997E+00,
			0.9275932609E+01, 0.4164915298E+00,
			0.3820341298E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.5218776196E+02, -0.1325278809E-01, 0.3759696623E-02,
			0.1032932006E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.3210344977E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.1235135428E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.5364201581E+00,  0.5951172526E+00, 0.4258595477E+00,
			0.2458806060E+00,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Flourine, BRACES(
		STO_6G_S(
			0.1728626574E+04, 0.9163596281E-02,
			0.3169417900E+03, 0.4936149294E-01,
			0.8866889139E+02, 0.1685383049E+00,
			0.3046015731E+02, 0.3705627997E+00,
			0.1182857044E+02, 0.4164915298E+00,
			0.4871658522E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.6703228091E+02, -0.1325278809E-01, 0.3759696623E-02,
			0.1326743777E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.4123509771E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.1586462839E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.6890018919E+00,  0.5951172526E+00, 0.4258595477E+00,
			0.3158199784E+00,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Neon, BRACES(
		STO_6G_S(
			0.2146955475E+04, 0.9163596281E-02,
			0.3936419362E+03, 0.4936149294E-01,
			0.1101268283E+03, 0.1685383049E+00,
			0.3783153777E+02, 0.3705627997E+00,
			0.1469109318E+02, 0.4164915298E+00,
			0.6050603466E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.8550442919E+02, -0.1325278809E-01, 0.3759696623E-02,
			0.1692355799E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.5259829210E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.2023645885E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.8788707870E+00,  0.5951172526E+00, 0.4258595477E+00,
			0.4028507849E+00,  0.2407061763E+00, 0.1017082955E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Sodium, BRACES(
		STO_6G_S(
			0.2600756771E+04, 0.9163596281E-02,
			0.4768459071E+03, 0.4936149294E-01,
			0.1334043010E+03, 0.1685383049E+00,
			0.4582797788E+02, 0.3705627997E+00,
			0.1779634487E+02, 0.4164915298E+00,
			0.7329517596E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.1248424044E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.2470956992E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.7679715913E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.2954663523E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.1283212382E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.5881901217E+00,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.9433006047E+01, -0.7943126362E-02, -0.7139358907E-02,
			0.2526243756E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.9473682506E+00, -0.1785026925E+00,  0.7621621429E-01,
			0.4240594435E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.2098454079E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1081470943E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Magnesium, BRACES(
		STO_6G_S(
			0.3103386324E+04, 0.9163596281E-02,
			0.5690024854E+03, 0.4936149294E-01,
			0.1591863906E+03, 0.1685383049E+00,
			0.5468482151E+02, 0.3705627997E+00,
			0.2123571643E+02, 0.4164915298E+00,
			0.8746040738E+01, 0.1303340841E+00
		),
		STO_6G_SP(
			0.1567952315E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.3103386828E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.9645303131E+01, -0.3378537151E-01, 0.1738967435E+00,
			0.3710895775E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.1611644564E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.7387346208E+00,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.8901677544E+01, -0.7943126362E-02, -0.7139358907E-02,
			0.2383949209E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.8940062839E+00, -0.1785026925E+00,  0.7621621429E-01,
			0.4001736462E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.1980255441E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1020555436E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Aluminum, BRACES(
		STO_6G_S(
			0.3644586388E+04, 0.9163596281E-02,
			0.6682309248E+03, 0.4936149294E-01,
			0.1869469321E+03, 0.1685383049E+00,
			0.6422131675E+02, 0.3705627997E+00,
			0.2493901660E+02, 0.4164915298E+00,
			0.1027126426E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.1959641441E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.3878641831E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.1205478990E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.4637912184E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.2014248422E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.9232774259E+00,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.8901677544E+01, -0.7943126362E-02, -0.7139358907E-02,
			0.2383949209E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.8940062839E+00, -0.1785026925E+00,  0.7621621429E-01,
			0.4001736462E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.1980255441E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1020555436E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Silicon, BRACES(
		STO_6G_S(
			0.4229261737E+04, 0.9163596281E-02,
			0.7754305100E+03, 0.4936149294E-01,
			0.2169375129E+03, 0.1685383049E+00,
			0.7452389069E+02, 0.3705627997E+00,
			0.2893980755E+02, 0.4164915298E+00,
			0.1191901091E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.2404904849E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.4759934318E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.1479384038E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.5691723632E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.2471919452E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.1133061545E+01,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.9433006047E+01, -0.7943126362E-02, -0.7139358907E-02,
			0.2526243756E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.9473682506E+00, -0.1785026925E+00,  0.7621621429E-01,
			0.4240594435E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.2098454079E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1081470943E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Phosphorous, BRACES(
		STO_6G_S(
			0.4857412371E+04, 0.9163596281E-02,
			0.8906012410E+03, 0.4936149294E-01,
			0.2491581331E+03, 0.1685383049E+00,
			0.8559254335E+02, 0.3705627997E+00,
			0.3323808927E+02, 0.4164915298E+00,
			0.1368928069E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.2906649590E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.5753018103E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.1788033738E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.6879210280E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.2987645712E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.1369456623E+01,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.1111939652E+02, -0.7943126362E-02, -0.7139358907E-02,
			0.2977874272E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.1116734493E+01, -0.1785026925E+00,  0.7621621429E-01,
			0.4998708868E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.2473606277E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1274811462E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Sulfur, BRACES(
		STO_6G_S(
			0.5529038289E+04, 0.9163596281E-02,
			0.1013743118E+04, 0.4936149294E-01,
			0.2836087927E+03, 0.1685383049E+00,
			0.9742727471E+02, 0.3705627997E+00,
			0.3783386178E+02, 0.4164915298E+00,
			0.1558207360E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.3455896791E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.6840121655E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.2125904712E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.8179121699E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.3552198128E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.1628232301E+01,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.1294439442E+02, -0.7943126362E-02, -0.7139358907E-02,
			0.3466625105E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.1300021248E+01, -0.1785026925E+00,  0.7621621429E-01,
			0.5819134077E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.2879592903E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1484042983E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Chlorine, BRACES(
		STO_6G_S(
			0.6236545525E+04, 0.9163596281E-02,
			0.1143463795E+04, 0.4936149294E-01,
			0.3198999635E+03, 0.1685383049E+00,
			0.1098942714E+03, 0.3705627997E+00,
			0.4267516141E+02, 0.4164915298E+00,
			0.1757598814E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.4039729660E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.7995679269E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.2485051157E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.9560887526E+01,  0.2502417861E+00, 0.4180364347E+00,
			0.4152299968E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.1903302881E+01,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.1358352871E+02, -0.7943126362E-02, -0.7139358907E-02,
			0.3637791008E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.1364210281E+01, -0.1785026925E+00,  0.7621621429E-01,
			0.6106455986E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.3021773873E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1557318157E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Argon, BRACES(
		STO_6G_S(
			0.6994673814E+04, 0.9163596281E-02,
			0.1282465787E+04, 0.4936149294E-01,
			0.3587877117E+03, 0.1685383049E+00,
			0.1232532624E+03, 0.3705627997E+00,
			0.4786284856E+02, 0.4164915298E+00,
			0.1971256419E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.4682992148E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.9268863609E+02, -0.4699171014E-01, 0.3767936984E-01,
			0.2880755901E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.1108330631E+02,  0.2502417861E+00, 0.4180364347E+00,
			0.4813487481E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.2206373495E+01,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.1672190907E+02, -0.7943126362E-02, -0.7139358907E-02,
			0.4478277461E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.1679401631E+01, -0.1785026925E+00,  0.7621621429E-01,
			0.7517310408E+00,  0.1510635058E+00,  0.4145098597E+00,
			0.3719933828E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.1917125747E+00,  0.2760593123E+00,  0.1058816521E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Potassium, BRACES(
		STO_6G_S(
			0.8001321412E+04, 0.9163596281E-02,
			0.1467033522E+04, 0.4936149294E-01,
			0.4104231128E+03, 0.1685383049E+00,
			0.1409914163E+03, 0.3705627997E+00,
			0.5475109279E+02, 0.4164915298E+00,
			0.2254952356E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.5433465051E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.1075424534E+03, -0.4699171014E-01, 0.3767936984E-01,
			0.3342411435E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.1285946155E+02,  0.2502417861E+00, 0.4180364347E+00,
			0.5584872913E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.2559955878E+01,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.2329374963E+02, -0.7943126362E-02, -0.7139358907E-02,
			0.6238275397E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.2339419558E+01, -0.1785026925E+00,  0.7621621429E-01,
			0.1047167197E+01,  0.1510635058E+00,  0.4145098597E+00,
			0.5181896807E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.2670571103E+00,  0.2760593123E+00,  0.1058816521E+00
		),
		STO_6G_SP(
			0.2791996035E+01,  0.3775056180E-02, -0.7052075733E-02,
			0.8983681264E+00, -0.5585965266E-01, -0.5259505547E-01,
			0.3838418398E+00, -0.3192946152E+00, -0.3773450392E-01,
			0.1914081612E+00, -0.2764780132E-01,  0.3874773403E+00,
			0.1033137261E+00,  0.9049199432E+00,  0.5791672602E+00,
			0.5744847995E-01,  0.3406258162E+00,  0.1221817127E+00
		)
	)),
	STO_6G_ATOM(ATOM_TYPE::Calcium, BRACES(
		STO_6G_S(
			0.8857157042E+04, 0.9163596281E-02,
			0.1623950048E+04, 0.4936149294E-01,
			0.4543227021E+03, 0.1685383049E+00,
			0.1560721100E+03, 0.3705627997E+00,
			0.6060736746E+02, 0.4164915298E+00,
			0.2496146087E+02, 0.1303340841E+00
		),
		STO_6G_SP(
			0.6175690999E+03, -0.1325278809E-01, 0.3759696623E-02,
			0.1222330419E+03, -0.4699171014E-01, 0.3767936984E-01,
			0.3798993832E+02, -0.3378537151E-01, 0.1738967435E+00,
			0.1461609860E+02,  0.2502417861E+00, 0.4180364347E+00,
			0.6347781583E+01,  0.5951172526E+00, 0.4258595477E+00,
			0.2909652740E+01,  0.2407061763E+00, 0.1017082955E+00
		),
		STO_6G_SP(
			0.2790660509E+02, -0.7943126362E-02, -0.7139358907E-02,
			0.7473639527E+01, -0.7100264172E-01, -0.1829277070E-01,
			0.2802694233E+01, -0.1785026925E+00,  0.7621621429E-01,
			0.1254537458E+01,  0.1510635058E+00,  0.4145098597E+00,
			0.6208066547E+00,  0.7354914767E+00,  0.4889621471E+00,
			0.3199423636E+00,  0.2760593123E+00,  0.1058816521E+00
		),
		STO_6G_SP(
			0.2525343962E+01,  0.3775056180E-02, -0.7052075733E-02,
			0.8125686765E+00, -0.5585965266E-01, -0.5259505547E-01,
			0.3471826822E+00, -0.3192946152E+00, -0.3773450392E-01,
			0.1731275539E+00, -0.2764780132E-01,  0.3874773403E+00,
			0.9344665645E-01,  0.9049199432E+00,  0.5791672602E+00,
			0.5196181158E-01,  0.3406258162E+00,  0.1221817127E+00
		)
	)),
}
};
}