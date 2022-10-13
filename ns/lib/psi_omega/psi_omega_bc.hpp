#ifndef PSI_OMEGA_Bc_HPP
#define PSI_OMEGA_Bc_HPP

#include "common.hpp"

struct APsiOmegaBc{
	virtual ~APsiOmegaBc() = default;
};

struct PsiOmegaBc_Input: APsiOmegaBc{
	PsiOmegaBc_Input(std::function<double(Point)> psifun): APsiOmegaBc(), psifun(psifun) {};
	const std::function<double(Point)> psifun;
};

struct PsiOmegaBc_Wall: APsiOmegaBc{
	PsiOmegaBc_Wall(double psival): APsiOmegaBc(), psival(psival) {};
	const double psival;
};

struct PsiOmegaBc_Outflow: APsiOmegaBc{
	PsiOmegaBc_Outflow(): APsiOmegaBc(){};
};

struct PsiOmegaBc_Symmetry: APsiOmegaBc{
	PsiOmegaBc_Symmetry(double psival): APsiOmegaBc(), psival(psival) {};
	const double psival;
};


#endif
