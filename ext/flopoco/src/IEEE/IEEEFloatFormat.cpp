/*
 * Floating Point FMA for FloPoCo
 *
 * Author: Florent de Dinechin
 * Copyright 2010 ENS-Lyon, INRIA, CNRS, UCBL
 *
 * This file is part of the FloPoCo project developed by the Arenaire team at
 * Ecole Normale Superieure de Lyon
 */

#include "IEEEFloatFormat.h"

using namespace std;
namespace flopoco
{
	static vector<IEEEFloatFormat> formats = {
	    IEEEFloatFormat("binary16", 5, 10),
	    IEEEFloatFormat("binary32", 8, 23),
	    IEEEFloatFormat("binary64", 11, 52),
	    IEEEFloatFormat("binary128", 15, 112),
	    IEEEFloatFormat("binary256", 19, 236),
	};

	IEEEFloatFormat::IEEEFloatFormat(const char *name, int wE, int wF)
	    : name(name), wE(wE), wF(wF)
	{
	}
	const vector<IEEEFloatFormat> &IEEEFloatFormat::getStandardFormats()
	{
		return formats;
	}

} // namespace flopoco