/*
  IEEE-compatible floating-point numbers for FloPoCo

  Author: F. de Dinechin

  This file is part of the FloPoCo project
  developed by the Arenaire team at Ecole Normale Superieure de Lyon
  
  Initial software.
  Copyright © ENS-Lyon, INRIA, CNRS, UCBL,  
  2008-2010.
  All rights reserved.

  */

#include "MPFRSetExp.hpp"
#include "utils.hpp"


namespace flopoco{


	MPFRSetExp::MPFRSetExp(mpfr_exp_t emin, mpfr_exp_t emax) {
		orig_emin = mpfr_get_emin();
		mpfr_set_emin(emin);

		orig_emax = mpfr_get_emax();
		mpfr_set_emax(emax);
	}

	MPFRSetExp MPFRSetExp::setupIEEE(int wE, int wF) {
		// emin and emax are specified for a mantissa in (0.5, 1)
		// The formula should evaluate to -1073 for doubles, see MPFR doc;
		int emin = -(1<<(wE-1)) - wF + 3; // -1024 - 52 + 3 
		// The formula should evaluate mpfr_t mp to 1024 for doubles, see MPFR doc;
		int emax = (1<<(wE-1));
		return MPFRSetExp(emin, emax);
	}

	MPFRSetExp::~MPFRSetExp() {
		mpfr_set_emin(orig_emin);
		mpfr_set_emax(orig_emax);
	}
	
}