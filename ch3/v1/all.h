#ifndef ALL_H
#define ALL_H

using type_calc = double;  //type used for calculations double, float, or something with better precision

namespace Const
{
	const double eps_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double q_e = 1.602176565e-19;		// C, electron charge
	const double amu = 1.660538921e-27;		// kg, atomic mass unit
	const double m_e = 9.10938215e-31;		// kg, electron mass
	const double k = 1.380648e-23;			// J/K, Boltzmann constant
	const double pi = 3.141592653;			// pi
	const double eV_to_K = q_e/k;				// 1eV in K ~ 11604
}

#endif