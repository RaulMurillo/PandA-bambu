
#include <iostream>
#include <sstream>

#include "gmp.h"
#include "mpfr.h"

#include "FixAtan2ByRecipMultAtan.hpp"
#include "IntMult/IntMultiplier.hpp"
#include "IntAddSubCmp/IntAdder.hpp"
#include "ConstMult/FixRealKCM.hpp"
#include "FixFunctions/FixFunctionByPiecewisePoly.hpp"
#include "FixFunctions/FixFunctionByMultipartiteTable.hpp"


using namespace std;

/* TODO Debugging:
There are still a few last-bit errors with the current setup
One is unavoidable, it is atan2 (0,0)
   TODO add a don't care to the test framework?
*/
namespace flopoco{
	// an an option for outputting the norm of the vector as well (scaled or not)



	FixAtan2ByRecipMultAtan::FixAtan2ByRecipMultAtan(OperatorPtr parentOp, Target* target_, int wIn_, int wOut_, int degree) :
		FixAtan2(parentOp, target_, wIn_, wOut_)
	{
		//int stage;
		srcFileName="FixAtan2ByRecipMultAtan";
		setCopyrightString ( "Matei Istoan, Florent de Dinechin (2012-...)" );
		useNumericStd_Unsigned();

#define MULTDEBUG 1
#if MULTDEBUG
		// The following blob is to remove when the multiplier version works
		// along with a setPlainVHDL near the end
		bool wasPlainVHDL=true;
		if(!getTarget()->plainVHDL()) {
			REPORT(0, "setting plainVHDL to true while we debug it");
			getTarget()->setPlainVHDL(true);
			wasPlainVHDL=false;
		}
#endif
		REPORT(DEBUG, "  degree=" << degree);

		ostringstream name;
		name << "FixAtan2ByRecipMultAtan_" << wIn_ << "_" << wOut_ << "_uid" << getNewUId();
		setNameWithFreqAndUID( name.str() );

		///////////// VHDL GENERATION

		// Range reduction code is shared, defined in FixAtan2
		buildQuadrantRangeReduction();
		buildScalingRangeReduction();

		//int sizeZ=wOut-2; // w-2 because two bits come from arg red

		//FixFunctionByPiecewisePoly* recipTable;
		Operator* recipTable;
		//FixFunctionByPiecewisePoly* atanTable;
		Operator* atanTable;
		int msbRecip, lsbRecip, msbProduct, lsbProduct, msbAtan, lsbAtan;
		msbAtan = -2; // bits 0 and -1 come from the range reduction
		lsbAtan = -wOut+1;
		msbRecip = 0; // 2/(1+x) in 0..1 for x in 0..1
		msbProduct = -1 ; // y/x between	0 and 1 but the faithful product may overflow a bit.
		if(degree==0) { // both tables are correctly rounded
			lsbRecip = -wOut+1; // see error analysis in the Arith2015 paper. It says we should have -w, but exhaustive test show that -w+1 work :)
			lsbProduct = -wOut+1; // It says we should have -w, but exhaustive test show that -w+1 work :)
		}
		else{ // both tables are faithful
			lsbRecip = -wOut; // see error analysis in the paper
			lsbProduct = -wOut;
		}
		// recip table computes 2/(1+x) once we have killed the MSB of XRS, which is always 1.
		vhdl << tab << declare("XRm1", wIn-2) << " <= XRS" << range(wIn-3,0)  << "; -- removing the MSB which is constantly 1" << endl;
		ostringstream invfun;
		invfun << "2/(1+x)-1b"<<lsbRecip;


#if 0 // Multipartite does not always work here and tends to die without grace, mostly a bug of the exploration in the implementation (for wIn < wOut)
		// TODO ! There is still a potential
		if(degree==1) {
			newInstance("FixFunctionByMultipartiteTable",
									join("reciprocal_uid", getNewUId()),
									"f="+invfun.str() + " signedIn=0 compressTIV=true lsbIn="+to_string(-wOut+2) + " lsbOut="+to_string(lsbRecip),
									"X=>XRm1",
									"Y=>R0");
				}
		else { }
#endif
			
		newInstance("FixFunctionByPiecewisePoly",
								join("reciprocal_uid", getNewUId()),
								"f="+invfun.str() + " signedIn=0 compressTIV=true lsbIn="+to_string(-wOut+2) + " lsbOut="+to_string(lsbRecip) + " d="+to_string(degree),
								"X=>XRm1",
								"Y=>R0");
		

		vhdl << tab << declareFixPoint("R", false, msbRecip, lsbRecip) << " <= unsigned(R0" << range(msbRecip-lsbRecip  , 0) << "); -- removing the sign  bit" << endl;
		vhdl << tab << declareFixPoint("YRU", false, -1, -wOut+1) << " <= unsigned(YRS);" << endl;

		if(getTarget()->plainVHDL()) { // generate a "*"
			vhdl << tab << declareFixPoint(getTarget()->DSPMultiplierDelay(), "P", false, msbRecip -1 +1, lsbRecip-wOut+1) << " <= R*YRU;" << endl;
			resizeFixPoint("PtruncU", "P", msbProduct, lsbProduct);
			vhdl << tab << declare("P_slv", msbProduct-lsbProduct+1)  << " <=  std_logic_vector(PtruncU);" << endl;
		}
		else{ // generate an IntMultiplier
			#if 0
			manageCriticalPath(getTarget()->DSPMultiplierDelay()); // This should be replaced with a method of IntMultiplier or something
			IntMultiplier::newComponentAndInstance(this,
																						 "divMult",     // instance name
																						 "R",  // x
																						 "YRU", // y
																						 "P",       // p
																						 msbProduct, lsbProduct
			);
			#else
			// TODO an instance
			#endif
		}


		string atanfun = "atan(x)/pi";
		
#if 0 // same comment as above, but now it dies because wIn>wOut: it seem we only ever tested wIn=wOut
		if(degree==1) {
			newInstance("FixFunctionByMultipartiteTable",
									join("atan_uid", getNewUId()),
									"f="+ atanfun+ " signedIn=0 compressTIV=true lsbIn="+to_string(lsbProduct) + " lsbOut="+to_string(lsbAtan),
									"X=>P_slv",
									"Y=>atanTableOut");
		}
		else {}
#endif
		newInstance("FixFunctionByPiecewisePoly",
								join("atan_uid", getNewUId()),
								"f="+ atanfun+ " signedIn=0 compressTIV=true lsbIn="+to_string(lsbProduct) + " lsbOut="+to_string(lsbAtan) + " d="+to_string(degree),
								"X=>P_slv",
								"Y=>atanTableOut");
	
		vhdl << tab << declare("finalZ", wOut) << " <= \"00\" & atanTableOut;" << endl;
		
		// Reconstruction code is shared, defined in FixAtan2
		buildQuadrantReconstruction();
		
#if MULTDEBUG
		getTarget()->setPlainVHDL(wasPlainVHDL);
#endif
	};


	FixAtan2ByRecipMultAtan::~FixAtan2ByRecipMultAtan(){
	};






}

