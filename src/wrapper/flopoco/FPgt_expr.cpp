/*
 *
 *                   _/_/_/    _/_/   _/    _/ _/_/_/    _/_/
 *                  _/   _/ _/    _/ _/_/  _/ _/   _/ _/    _/
 *                 _/_/_/  _/_/_/_/ _/  _/_/ _/   _/ _/_/_/_/
 *                _/      _/    _/ _/    _/ _/   _/ _/    _/
 *               _/      _/    _/ _/    _/ _/_/_/  _/    _/
 *
 *             ***********************************************
 *                              PandA Project
 *                     URL: http://panda.dei.polimi.it
 *                       Politecnico di Milano - DEIB
 *                        System Architectures Group
 *             ***********************************************
 *              Copyright (C) 2004-2020 Politecnico di Milano
 *
 *   This file is part of the PandA framework.
 *
 *   The PandA framework is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
/**
 * @file FPgt_expr.cpp
 * @brief FPgt_expr module for flopoco.
 *
 * @author Fabrizio Ferrandi <fabrizio.ferrandi@polimi.it>
 * $Date$
 * Last modified by $Author$
 *
 */

/// Autoheader include
#include "config_SKIP_WARNING_SECTIONS.hpp"

#if SKIP_WARNING_SECTIONS
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <cmath>
#include <cstring>
#include <iosfwd>
#include <sstream>
#include <vector>

#include <cstddef>
#include <gmp.h>

#include "utils.hpp"
#include <gmpxx.h>

#include "FPAddSub/FPAddSinglePath.hpp"
#include "FPgt_expr.hpp"

#include "custom_map.hpp"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <list>
#include <locale>
#include <sstream>
#include <string>
#include <vector>

#include <cstdio>
#include <mpfr.h>

#include "flopoco_wrapper.hpp"

using namespace std;

namespace flopoco
{

#define DEBUGVHDL 0

   FPgt_expr::FPgt_expr(Operator* parentOp, Target* _target, int wE, int wF) : Operator(parentOp, _target)
   {
      ostringstream name;

      name << "FPgt_expr_" << wE << "_" << wF;
      setName(name.str());

      setCopyrightString("Fabrizio Ferrandi (2011-2018)");

      /* Set up the IO signals */

      addFPInput("X", wE, wF);
      addFPInput("Y", wE, wF);
      addOutput("R");

      /*	VHDL code description	*/
      vhdl << tab << declare(getTarget()->logicDelay(1), "nY", wE + wF + 3) << "  <= Y" << range(wE + wF + 2, wE + wF + 1) << " & not(Y" << of(wE + wF) << ") & Y" << range(wE + wF - 1, 0) << ";" << endl;

      ostringstream paramR, inmapR, outmapR;
      paramR << "wE=" << wE;
      paramR << " wF=" << wF;

      inmapR << "X=>X, Y=>nY";
      outmapR << "R=>valueDiff";

      newInstance("FPAddSinglePath", "value_difference", paramR.str(), inmapR.str(), outmapR.str());

      vhdl << tab << declare(getTarget()->logicDelay(2), "R0", 1, false) << " <= '1' when (valueDiff" << of(wE + wF) << "='0') and (valueDiff" << range(wE + wF + 2, wE + wF + 1) << " /= \"00\") else '0';" << endl;
      vhdl << tab << "R <= R0;" << endl;
   }

   FPgt_expr::~FPgt_expr() = default;

   void FPgt_expr::emulate(TestCase*)
   {
      // TODO
   }

   void FPgt_expr::buildStandardTestCases(TestCaseList*)
   {
   }

} // namespace flopoco
