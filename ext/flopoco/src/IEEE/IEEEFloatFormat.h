#ifndef IEEE_FLOAT_FORMAT_HPP
#define IEEE_FLOAT_FORMAT_HPP
#include <vector>

namespace flopoco
{
	/**
	 * Define IEEE754 floating point standard formats
	 */
	class IEEEFloatFormat
	{
	public:
		static const std::vector<IEEEFloatFormat> &getStandardFormats();
		const char *name;
		int wE; /** Exponent bits */
		int wF; /** Fraction bits */

		IEEEFloatFormat(const char *name, int wE, int wF);
	};
} // namespace flopoco

#endif