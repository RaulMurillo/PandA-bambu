#ifndef __MPFRSETEXP_HPP
#define __MPFRSETEXP_HPP

#include <gmpxx.h>
#include <mpfr.h>


namespace flopoco{

	/**
	 * A helper class to set and restore mpfr emin/emax.
	 */
	class MPFRSetExp
	{
	public:
		MPFRSetExp(mpfr_exp_t emin, mpfr_exp_t emax);
		~MPFRSetExp();

		/** Setup mpfr emin/emax for IEEE floating point format */
		static MPFRSetExp setupIEEE(int wE, int wF);

	private:
		mpfr_exp_t orig_emin;
		mpfr_exp_t orig_emax;

	};

}

#endif

