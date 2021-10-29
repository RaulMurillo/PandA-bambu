#ifndef __FPLOG_HPP
#define __FPLOG_HPP
#include <vector>
#include <sstream>

#include "Operator.hpp"

namespace flopoco{

#define MAXRRSTAGES 2000 // 4000 bits of accuracy should be enough for anybody

	class FPLog : public Operator
	{
	public:
		FPLog(OperatorPtr parentOp, Target* target, int wE, int wF);
		~FPLog();

		//		Overloading the virtual functions of Operator
		void emulate(TestCase * tc);
		void buildStandardTestCases(TestCaseList* tcl);
		static TestList unitTest(int index);
		/**Overloading the function of Operator with a function that tests only positive FP numbers (full range)*/
		TestCase* buildRandomTestCase(int i);
		// User-interface stuff
		/** Factory method */
		static OperatorPtr parseArguments(OperatorPtr parentOp, Target *target , vector<string> &args);
		static void registerFactory();
	protected:
		int wE, wF;

	};
	
}
#endif
