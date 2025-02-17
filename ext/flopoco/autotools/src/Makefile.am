# Makefile for flopoco/src
SUBDIRS = 
DIST_SUBDIRS =

EXTRA_DIST = VHDLLexerRule.l Factories/MakeFactory.sh FPExpressions/ExpressionScanner.l FPExpressions/ExpressionParser.y \
		$(top_srcdir)/src/FPExpressions/ExpressionParser.h $(top_srcdir)/src/FPExpressions/FlopocoScanner.h
ACLOCAL_AMFLAGS = -I ../m4

BUILT_SOURCES = VHDLLexer.cpp Factories.cpp ExpressionParser.cpp ExpressionScanner.cpp
                
CLEANFILES = $(BUILT_SOURCES) lex.yy.c *.output

noinst_LTLIBRARIES = libflopoco.la

noinst_HEADERS = \
		Operator.hpp \
		UserInterface.hpp \
		Signal.hpp \
		Target.hpp \
		utils.hpp \
		FlopocoStream.hpp \
		TargetModel.hpp \
		TutorialOperator.hpp \
		Instance.hpp \
		ShiftReg.hpp \
		Complex/FPComplexMultiplier.hpp \
		Complex/OldFixedComplexMultiplier.hpp \
		Complex/FixComplexMult.hpp \
		Complex/FixFFTFullyPA.hpp \
		Complex/IntTwiddleMultiplier.hpp \
		Complex/FixFFT.hpp \
		Complex/FixComplexR2Butterfly.hpp \
		Complex/IntFFTLevelDIT2.hpp \
		Complex/FPComplexAdder.hpp \
		Complex/FixComplexKCM.hpp \
		Complex/IntTwiddleMultiplierAlternative.hpp \
		Complex/FixComplexAdd.hpp \
		Complex/IntFFT.hpp \
		Complex/IntFFTButterfly.hpp \
		Complex/FixComplexAdder.hpp \
		FixFunctions/VaryingPiecewisePolyApprox.hpp \
		FixFunctions/FixFunctionByMultipartiteTable.hpp \
		FixFunctions/UniformPiecewisePolyApprox.hpp \
		FixFunctions/BasicPolyApprox.hpp \
		FixFunctions/HOTBM.hpp \
		FixFunctions/HOTBM/TermPowMult.hh \
		FixFunctions/HOTBM/PWPolynomial.hh \
		FixFunctions/HOTBM/Component.hh \
		FixFunctions/HOTBM/PowerAdHoc.hh \
		FixFunctions/HOTBM/Util.hh \
		FixFunctions/HOTBM/TermROM.hh \
		FixFunctions/HOTBM/Minimax.hh \
		FixFunctions/HOTBM/Exhaustive.hh \
		FixFunctions/HOTBM/PowerROM.hh \
		FixFunctions/HOTBM/Polynomial.hh \
		FixFunctions/HOTBM/MPPolynomial.hh \
		FixFunctions/HOTBM/HOTBMInstance.hh \
		FixFunctions/HOTBM/Term.hh \
		FixFunctions/HOTBM/Param.hh \
		FixFunctions/HOTBM/Power.hh \
		FixFunctions/FixFunctionBySimplePoly.hpp \
		FixFunctions/FixFunctionByVaryingPiecewisePoly.hpp \
		FixFunctions/Attic/PolyCoeffTable.hpp \
		FixFunctions/Attic/FunctionEvaluator.hpp \
		FixFunctions/Attic/PiecewiseFunction.hpp \
		FixFunctions/Attic/PolynomialEvaluator.hpp \
		FixFunctions/Attic/BipartiteTable.hpp \
		FixFunctions/Attic/Function.hpp \
		FixFunctions/FixPolyEval.hpp \
		FixFunctions/FixFunctionByTable.hpp \
		FixFunctions/FixFunction.hpp \
		FixFunctions/FixHornerEvaluator.hpp \
		FixFunctions/FixFunctionByPiecewisePoly.hpp \
		FixFunctions/Multipartite.hpp \
		ConstMult/CRFPConstMult.hpp \
		ConstMult/adder_cost.hpp \
		ConstMult/IntConstMultShiftAdd.hpp \
		ConstMult/IntIntKCM.hpp \
		ConstMult/error_comp_graph.hpp \
		ConstMult/FPConstMult.hpp \
		ConstMult/IntConstMultShiftAddRPAG.hpp \
		ConstMult/FixRealKCM.hpp \
		ConstMult/tscm_solutions.hpp \
		ConstMult/WordLengthCalculator.hpp \
		ConstMult/FPConstDiv.hpp \
		ConstMult/ShiftAddDag.hpp \
		ConstMult/FixFixConstMult.hpp \
		ConstMult/FPRealKCM.hpp \
		ConstMult/IntConstMultShiftAddTypes.hpp \
		ConstMult/IntConstMultShiftAddOptTernary.hpp \
		ConstMult/ShiftAddOp.hpp \
		ConstMult/FixRealConstMult.hpp \
		ConstMult/KCMTable.hpp \
		ConstMult/IntConstDiv.hpp \
		ConstMult/IntConstMCM.hpp \
		ConstMult/FixRealShiftAdd.hpp \
		ConstMult/IntConstMultShiftAddOpt.hpp \
		ConstMult/IntConstMult.hpp \
		LexerContext.hpp \
		FixFilters/FixIIRShiftAdd.hpp \
		FixFilters/FixRootRaisedCosine.hpp \
		FixFilters/Attic/FixDCT.hpp \
		FixFilters/Attic/FixRCF.hpp \
		FixFilters/FixFIR.hpp \
		FixFilters/FixHalfSine.hpp \
		FixFilters/FixIIR.hpp \
		FixFilters/FixSOPC.hpp \
		FPDivSqrt/FPDiv.hpp \
		FPDivSqrt/FPSqrt.hpp \
		FPDivSqrt/FPSqrtPoly.hpp \
		FPDivSqrt/Attic/PolynomialTableCorrectRounded.hpp \
		FPDivSqrt/Attic/PolynomialTable.hpp \
		Trigs/FixSinPoly.hpp \
		Trigs/FixSinCos.hpp \
		Trigs/FixSinCosCORDIC.hpp \
		Trigs/FixAtan2.hpp \
		Trigs/FixSinOrCos.hpp \
		Trigs/Atan2Table.hpp \
		Trigs/FixAtan2ByCORDIC.hpp \
		Trigs/FixAtan2ByRecipMultAtan.hpp \
		Trigs/FixSinCosPoly.hpp \
		Trigs/Fix2DNorm.hpp \
		Trigs/FixAtan2ByBivariateApprox.hpp \
		Trigs/ConstDiv3ForSinPoly.hpp \
		OperatorContext.hpp \
		PrimitiveComponents/Intel/Intel_LCELL.hpp \
		PrimitiveComponents/GenericAddSub.hpp \
		PrimitiveComponents/GenericMult.hpp \
		PrimitiveComponents/BooleanEquation.hpp \
		PrimitiveComponents/GenericMux.hpp \
		PrimitiveComponents/Primitive.hpp \
		PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State_slice.hpp \
		PrimitiveComponents/Xilinx/Xilinx_MUXF7.hpp \
		PrimitiveComponents/Xilinx/Xilinx_LUT6.hpp \
		PrimitiveComponents/Xilinx/Xilinx_n2m_Decoder.hpp \
		PrimitiveComponents/Xilinx/XilinxGPC.hpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericMux_slice.hpp \
		PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State.hpp \
		PrimitiveComponents/Xilinx/Xilinx_Comparator.hpp \
		PrimitiveComponents/Xilinx/Xilinx_CFGLUT5.hpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericAddSub_slice.hpp \
		PrimitiveComponents/Xilinx/Xilinx_CARRY4.hpp \
		PrimitiveComponents/Xilinx/Xilinx_Primitive.hpp \
		PrimitiveComponents/Xilinx/Xilinx_LUT5.hpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericAddSub.hpp \
		PrimitiveComponents/Xilinx/Xilinx_MUXF8.hpp \
		PrimitiveComponents/Xilinx/XilinxFourToTwoCompressor.hpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericMux.hpp \
		PrimitiveComponents/Xilinx/Xilinx_FDCE.hpp \
		PrimitiveComponents/Xilinx/Xilinx_CFGLUTShadow.hpp \
		PrimitiveComponents/Xilinx/Xilinx_LUT_compute.h \
		PrimitiveComponents/GenericLut.hpp \
		FloPoCo.hpp \
		Conversions/InputIEEE.hpp \
		Conversions/OutputIEEE.hpp \
		Conversions/FP2Fix.hpp \
		Conversions/Fix2FP.hpp \
		Conversions/PIF2Posit.hpp \
		Conversions/Posit2PIF.hpp \
		Conversions/FP2Posit.hpp \
		Conversions/Posit2FP.hpp \
		Conversions/Posit2Int.hpp \
		Conversions/Posit2UInt.hpp \
		Conversions/Int2Posit.hpp \
		Conversions/UInt2Posit.hpp \
		Conversions/Posit2Quire.hpp \
		Conversions/Quire2Posit.hpp \
		Table.hpp \
		main.hpp \
		ShiftersEtc/Shifters.hpp \
		ShiftersEtc/Normalizer.hpp \
		ShiftersEtc/LZOC3.hpp \
		ShiftersEtc/LZOC.hpp \
		BitHeap/OptimalCompressionStrategy.hpp \
		BitHeap/Bit.hpp \
		BitHeap/MaxEfficiencyCompressionStrategy.hpp \
		BitHeap/BitHeap.hpp \
		BitHeap/WeightedBit.hpp \
		BitHeap/BitHeapPlotter.hpp \
		BitHeap/ParandehAfsharCompressionStrategy.hpp \
		BitHeap/BitHeapTest.hpp \
		BitHeap/FirstFittingCompressionStrategy.hpp \
		BitHeap/Compressor.hpp \
		BitHeap/Attic/BitHeap.hpp \
		BitHeap/Attic/Plotter.hpp \
		BitHeap/BitHeapSolution.hpp \
		BitHeap/CompressionStrategy.hpp \
		Targets/AllTargetsHeaders.hpp \
		Targets/Kintex7.hpp \
		Targets/Zynq7000.hpp \
		Targets/Virtex6.hpp \
		Targets/VirtexUltrascalePlus.hpp \
		Targets/DSP.hpp \
		Targets/StratixV.hpp \
		Targets/Old/StratixIV.hpp \
		Targets/Old/CycloneIII.hpp \
		Targets/Old/StratixII.hpp \
		Targets/Old/CycloneII.hpp \
		Targets/Old/StratixIII.hpp \
		Targets/Old/Virtex4.hpp \
		Targets/Old/Spartan3.hpp \
		Targets/Old/CycloneIV.hpp \
		Targets/Old/CycloneV.hpp \
		Targets/Old/Virtex5.hpp \
		AutoTest/AutoTest.hpp \
		random/moment_correction/systems/symmetric_quintic_system.hpp \
		random/moment_correction/systems/symmetric_heptic_system.hpp \
		random/moment_correction/systems/cubic_system.hpp \
		random/moment_correction/systems/quintic_system.hpp \
		random/moment_correction/systems/symmetric_nonic_system.hpp \
		random/moment_correction/systems/symmetric_cubic_system.hpp \
		random/moment_correction/correct_distribution.hpp \
		random/distributions/table_distribution.hpp \
		random/distributions/moment_conversions.hpp \
		random/distributions/distribution.hpp \
		random/distributions/make_table_approximation.hpp \
		random/distributions/gaussian_distribution.hpp \
		random/urng/LutSrRng.hpp \
		random/utils/minimise_gsl.hpp \
		random/utils/gsl_utils.hpp \
		random/utils/conversions.hpp \
		random/utils/find_roots_gsl.hpp \
		OperatorPipeline/structure_donnees/File.hpp \
		OperatorPipeline/structure_donnees/OPOperation.hpp \
		OperatorPipeline/structure_donnees/OPVariable.hpp \
		OperatorPipeline/structure_donnees/OPFunction.hpp \
		OperatorPipeline/structure_donnees/OPExpression.hpp \
		OperatorPipeline/OperatorPipeline.hpp \
		OperatorPipeline/ExpressionTreeData.h \
		IntAddSubCmp/IntDualAddSub.hpp \
		IntAddSubCmp/IntComparator.hpp \
		IntAddSubCmp/CarryGenerationCircuit.hpp \
		IntAddSubCmp/Attic/LongIntAdderCmpAddIncGen1.hpp \
		IntAddSubCmp/Attic/LongIntAdderAddAddMuxGen2.hpp \
		IntAddSubCmp/Attic/LongIntAdderMuxNetwork.hpp \
		IntAddSubCmp/Attic/LongIntAdderCmpCmpAddGen1.hpp \
		IntAddSubCmp/Attic/IntAdderAlternative.hpp \
		IntAddSubCmp/Attic/IntAdderSpecific.hpp \
		IntAddSubCmp/Attic/LongIntAdderCmpCmpAddGen2.hpp \
		IntAddSubCmp/Attic/LongIntAdderAddAddMuxGen1.hpp \
		IntAddSubCmp/Attic/LongIntAdderCmpAddIncGen2.hpp \
		IntAddSubCmp/Attic/IntAdderClassical.hpp \
		IntAddSubCmp/Attic/IntAdderShortLatency.hpp \
		IntAddSubCmp/IntMultiAdder.hpp \
		IntAddSubCmp/IntComparatorSpecific.hpp \
		IntAddSubCmp/IntAdder.hpp \
		Permuters/Permuter.hpp \
		FPExpressions/FlopocoScanner.h \
		FPExpressions/ExpressionParser.h \
		FPExpressions/ExpressionParserData.h \
		ExpLog/FPExpAttic/fragment.h \
		ExpLog/FPExpAttic/logfragment.h \
		ExpLog/FPExpAttic/gen_table.h \
		ExpLog/FPExpAttic/explore.h \
		ExpLog/FPExpAttic/math_lib.h \
		ExpLog/FPPow.hpp \
		ExpLog/FPExp.hpp \
		ExpLog/FPLog.hpp \
		ExpLog/FPLogIterative.hpp \
		FPComposite/FPSumOf3Squares.hpp \
		FPComposite/FPLargeAcc.hpp \
		FPComposite/LargeAccToFP.hpp \
		FPComposite/FPDotProduct.hpp \
		FPComposite/FP2DNorm.hpp \
		Attic/FPPipeline.hpp \
		IntMult/FixMultAdd.hpp \
		IntMult/BaseMultiplierDSPSuperTilesXilinx.hpp \
		IntMult/BaseMultiplier.hpp \
		IntMult/IntSquarer.hpp \
		IntMult/BaseMultiplierDSP.hpp \
		IntMult/IntMultiplierLUT.hpp \
		IntMult/IntKaratsuba.hpp \
		IntMult/BaseMultiplierLUT.hpp \
		IntMult/FormalBinaryProduct.hpp \
		IntMult/BaseMultiplierIrregularLUTXilinx.hpp \
		IntMult/Tiling.hpp \
		IntMult/BaseMultiplierCollection.hpp \
		IntMult/MultiplierTileCollection.hpp \
		IntMult/TilingStrategy.hpp \
		IntMult/TilingAndCompressionOptILP.hpp \
		IntMult/BaseMultiplierDSPKaratsuba.hpp \
		IntMult/LineCursor.hpp \
		IntMult/TilingStrategyGreedy.hpp \
		IntMult/Attic/FixMultAdd.hpp \
		IntMult/Attic/IntMultiplierOld.hpp \
		IntMult/IntMultiplier.hpp \
		IntMult/DSPBlock.hpp \
		IntMult/MultiplierBlock.hpp \
		IntMult/TilingStrategyXGreedy.hpp \
		IntMult/TilingStrategyBasicTiling.hpp \
		IntMult/NearestPointCursor.hpp \
		IntMult/GenericBinaryPolynomial.hpp \
		IntMult/TargetBaseMultSpecialisation.hpp \
		IntMult/IntKaratsubaRectangular.hpp \
		IntMult/TilingStrategyBeamSearch.hpp \
		IntMult/BaseMultiplierCategory.hpp \
		IntMult/BaseMultiplierXilinx2xk.hpp \
		IntMult/Field.hpp \
		IntMult/TilingStrategyOptimalILP.hpp \
		IntMult/IntPower.hpp \
		LNS/CotranTables.hpp \
		LNS/CotranHybrid.hpp \
		LNS/LNSMul.hpp \
		LNS/LNSDiv.hpp \
		LNS/LNSAdd.hpp \
		LNS/LogSinCos.hpp \
		LNS/LNSSqrt.hpp \
		LNS/LNSAddSub.hpp \
		LNS/Cotran.hpp \
		LNS/AtanPow.hpp \
		FPAddSub/FPAddSinglePath.hpp \
		FPAddSub/FPAdd3Input.hpp \
		FPAddSub/FPAddDualPath.hpp \
		FPAddSub/FPAdd.hpp \
		FPAddSub/FPAddSub.hpp \
		IEEE/IEEEAdd.hpp \
		IEEE/IEEEFMA.hpp \
		IEEE/IEEEFloatFormat.h \
		Posit/Posit2Posit.hpp \
		Posit/Add/PIFAdd.hpp \
		Posit/Add/PositAdd.hpp \
		Posit/Add/PositAdder.hpp \
		Posit/Add/PositAddSub.hpp \
		Posit/Exp/PIF2Fix.hpp \
		Posit/Exp/PositExp.hpp \
		Posit/Fun/PositFunction.hpp \
		Posit/Fun/PositFunctionByTable.hpp \
		Posit/Decoder/PositDecoder.hpp \
		Posit/Decoder/PositFastDecoder.hpp \
		Posit/Encoder/PositEncoder.hpp \
		Posit/Encoder/PositFastEncoder.hpp \
		Posit/Mult/PositMult.hpp \
		Posit/MAC/PositMAC.hpp \
		Posit/Div/PositDiv.hpp \
		Posit/Sqrt/PositSqrt.hpp \
		Posit/ApproxDiv/PositApproxDiv.hpp \
		Posit/LogApproxMult/PositLAM.hpp \
		Posit/Cmp/PositComparator.hpp \
		Posit/PositAssign.hpp \
		DualTable.hpp \
		TestBenches/Wrapper.hpp \
		TestBenches/IEEENumber.hpp \
		TestBenches/TestBench.hpp \
		TestBenches/TestCase.hpp \
		TestBenches/PositNumber.hpp \
		TestBenches/FPNumber.hpp \
		TestBenches/MPFRSetExp.hpp \
		FPMultSquare/FPMult.hpp \
		FPMultSquare/FPSquare.hpp \
		FPMultSquare/FPMultKaratsuba.hpp \
		Tools/ResourceEstimationHelper.hpp \
		Tools/FloorplanningHelper.hpp \
		Tools/Plane.hpp \
		Tools/Point.hpp \
		FixFunctions/FixDiv.hpp \
		FixConstant.hpp

libflopoco_la_SOURCES = \
		Operator.cpp \
		UserInterface.cpp \
		Signal.cpp \
		Target.cpp \
		utils.cpp \
		FlopocoStream.cpp \
		VHDLLexer.cpp \
		Factories.cpp \
		Instance.cpp \
		Tools/ResourceEstimationHelper.cpp \
		Tools/FloorplanningHelper.cpp \
		Targets/DSP.cpp \
		Targets/Virtex6.cpp \
		Targets/Zynq7000.cpp \
		Targets/Kintex7.cpp \
		Targets/VirtexUltrascalePlus.cpp \
		Targets/StratixV.cpp \
		AutoTest/AutoTest.cpp \
		TestBenches/TestCase.cpp \
		TestBenches/FPNumber.cpp \
		TestBenches/IEEENumber.cpp \
		TestBenches/PositNumber.cpp \
		TestBenches/Wrapper.cpp \
		TestBenches/TestBench.cpp \
		TestBenches/MPFRSetExp.cpp \
		BitHeap/Bit.cpp \
		BitHeap/Compressor.cpp \
		BitHeap/BitHeap.cpp \
		BitHeap/CompressionStrategy.cpp \
		BitHeap/BitHeapPlotter.cpp \
		BitHeap/BitHeapSolution.cpp \
		BitHeap/BitHeapTest.cpp \
		BitHeap/FirstFittingCompressionStrategy.cpp \
		BitHeap/ParandehAfsharCompressionStrategy.cpp \
		BitHeap/MaxEfficiencyCompressionStrategy.cpp \
		BitHeap/OptimalCompressionStrategy.cpp \
		TutorialOperator.cpp \
		ShiftersEtc/LZOC.cpp \
		ShiftersEtc/LZOC3.cpp \
		ShiftersEtc/Normalizer.cpp \
		ShiftersEtc/Shifters.cpp \
		ShiftReg.cpp \
		FixFilters/FixSOPC.cpp \
		FixFilters/FixFIR.cpp \
		FixFilters/FixHalfSine.cpp \
		FixFilters/FixRootRaisedCosine.cpp \
		FixFilters/FixIIR.cpp \
		FixFilters/FixIIRShiftAdd.cpp \
		IntAddSubCmp/IntAdder.cpp \
		IntAddSubCmp/IntDualAddSub.cpp \
		IntAddSubCmp/IntMultiAdder.cpp \
		ConstMult/ShiftAddOp.cpp \
		ConstMult/ShiftAddDag.cpp \
		ConstMult/FixRealKCM.cpp \
		ConstMult/FixRealShiftAdd.cpp \
		ConstMult/FixRealConstMult.cpp \
		ConstMult/FixFixConstMult.cpp \
		ConstMult/IntConstMult.cpp \
		ConstMult/FPConstMult.cpp \
		ConstMult/IntConstDiv.cpp \
		ConstMult/IntConstMultShiftAdd.cpp \
		ConstMult/IntConstMultShiftAddTypes.cpp \
		ConstMult/IntConstMultShiftAddOptTernary.cpp \
		ConstMult/IntConstMultShiftAddOpt.cpp \
		ConstMult/IntConstMultShiftAddRPAG.cpp \
		ConstMult/adder_cost.cpp \
		ConstMult/error_comp_graph.cpp \
		ConstMult/WordLengthCalculator.cpp \
		IntMult/IntMultiplier.cpp \
		IntMult/FixMultAdd.cpp \
		IntMult/TilingStrategy.cpp \
		IntMult/TilingStrategyOptimalILP.cpp \
		IntMult/TilingStrategyBasicTiling.cpp \
		IntMult/TilingStrategyGreedy.cpp \
		IntMult/TilingStrategyXGreedy.cpp \
		IntMult/TilingStrategyBeamSearch.cpp \
		IntMult/Field.cpp \
		IntMult/NearestPointCursor.cpp \
		IntMult/MultiplierBlock.cpp \
		IntMult/DSPBlock.cpp \
		IntMult/BaseMultiplier.cpp \
		IntMult/BaseMultiplierCategory.cpp \
		IntMult/BaseMultiplierCollection.cpp \
		IntMult/BaseMultiplierDSP.cpp \
		IntMult/BaseMultiplierDSPKaratsuba.cpp \
		IntMult/BaseMultiplierDSPSuperTilesXilinx.cpp \
		IntMult/BaseMultiplierLUT.cpp \
		IntMult/IntKaratsubaRectangular.cpp \
		IntMult/IntMultiplierLUT.cpp \
		IntMult/TargetBaseMultSpecialisation.cpp \
		IntMult/BaseMultiplierXilinx2xk.cpp \
		IntMult/BaseMultiplierIrregularLUTXilinx.cpp \
		IntMult/MultiplierTileCollection.cpp \
		IntMult/BaseMultiplierDSPKaratsuba.cpp \
		IntMult/TilingAndCompressionOptILP.cpp \
		IntMult/IntSquarer.cpp \
		Complex/FixComplexKCM.cpp \
		Complex/FixComplexAdder.cpp \
		Complex/FixComplexR2Butterfly.cpp \
		Complex/FixFFTFullyPA.cpp \
		Conversions/FP2Fix.cpp \
		Conversions/Fix2FP.cpp \
		Conversions/InputIEEE.cpp \
		Conversions/OutputIEEE.cpp \
		FPAddSub/FPAdd.cpp \
		FPAddSub/FPAddDualPath.cpp \
		FPAddSub/FPAddSinglePath.cpp \
		FPMultSquare/FPMult.cpp \
		FPDivSqrt/FPDiv.cpp \
		FPDivSqrt/FPSqrt.cpp \
		ExpLog/FPExp.cpp \
		ExpLog/FPLog.cpp \
		ExpLog/FPLogIterative.cpp \
		ExpLog/FPPow.cpp \
		IEEE/IEEEAdd.cpp \
		IEEE/IEEEFMA.cpp \
		IEEE/IEEEFloatFormat.cpp \
		Table.cpp \
		DualTable.cpp \
		FixConstant.cpp \
		FixFunctions/FixFunction.cpp \
		FixFunctions/FixFunctionByTable.cpp \
		FixFunctions/BasicPolyApprox.cpp \
		FixFunctions/FixPolyEval.cpp \
		FixFunctions/FixHornerEvaluator.cpp \
		FixFunctions/UniformPiecewisePolyApprox.cpp \
		FixFunctions/FixFunctionBySimplePoly.cpp \
		FixFunctions/FixFunctionByPiecewisePoly.cpp \
		FixFunctions/Multipartite.cpp \
		FixFunctions/FixFunctionByMultipartiteTable.cpp \
		Trigs/ConstDiv3ForSinPoly.cpp \
		Trigs/FixSinPoly.cpp \
		Trigs/FixSinCos.cpp \
		Trigs/FixSinCosPoly.cpp \
		Trigs/FixSinCosCORDIC.cpp \
		Trigs/FixAtan2.cpp \
		Trigs/FixAtan2ByCORDIC.cpp \
		Trigs/Fix2DNorm.cpp \
		Trigs/FixAtan2ByRecipMultAtan.cpp \
		PrimitiveComponents/BooleanEquation.cpp \
		PrimitiveComponents/GenericLut.cpp \
		PrimitiveComponents/Primitive.cpp \
		PrimitiveComponents/GenericAddSub.cpp \
		PrimitiveComponents/GenericMux.cpp \
		PrimitiveComponents/Xilinx/Xilinx_Primitive.cpp \
		PrimitiveComponents/Xilinx/Xilinx_LUT_compute.cpp \
		PrimitiveComponents/Xilinx/Xilinx_CARRY4.cpp \
		PrimitiveComponents/Xilinx/Xilinx_LUT5.cpp \
		PrimitiveComponents/Xilinx/Xilinx_LUT6.cpp \
		PrimitiveComponents/Xilinx/Xilinx_MUXF7.cpp \
		PrimitiveComponents/Xilinx/Xilinx_MUXF8.cpp \
		PrimitiveComponents/Xilinx/Xilinx_FDCE.cpp \
		PrimitiveComponents/Xilinx/Xilinx_n2m_Decoder.cpp \
		PrimitiveComponents/Xilinx/Xilinx_Comparator.cpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericAddSub.cpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericAddSub_slice.cpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericMux.cpp \
		PrimitiveComponents/Xilinx/Xilinx_GenericMux_slice.cpp \
		PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State.cpp \
		PrimitiveComponents/Xilinx/Xilinx_TernaryAdd_2State_slice.cpp \
		PrimitiveComponents/Xilinx/XilinxGPC.cpp \
		PrimitiveComponents/Xilinx/XilinxFourToTwoCompressor.cpp \
		PrimitiveComponents/Intel/Intel_LCELL.cpp \
		Tools/Plane.cpp \
		Tools/Point.cpp \
		Complex/FixFFT.cpp \
		Conversions/PIF2Posit.cpp \
		Conversions/Posit2PIF.cpp \
		Conversions/FP2Posit.cpp \
		Conversions/Posit2FP.cpp \
		Conversions/Posit2Int.cpp \
		Conversions/Posit2UInt.cpp \
		Conversions/Int2Posit.cpp \
		Conversions/UInt2Posit.cpp \
		Conversions/Posit2Quire.cpp \
		Conversions/Quire2Posit.cpp \
		Posit/Add/PIFAdd.cpp \
		Posit/Add/PositAdd.cpp \
		Posit/Add/PositAdder.cpp \
		Posit/Add/PositAddSub.cpp \
		Posit/Exp/PIF2Fix.cpp \
		Posit/Exp/PositExp.cpp \
		Posit/Fun/PositFunction.cpp \
		Posit/Fun/PositFunctionByTable.cpp \
		Posit/Posit2Posit.cpp \
		Posit/Decoder/PositDecoder.cpp \
		Posit/Decoder/PositFastDecoder.cpp \
		Posit/Encoder/PositEncoder.cpp \
		Posit/Encoder/PositFastEncoder.cpp \
		Posit/Mult/PositMult.cpp \
		Posit/MAC/PositMAC.cpp \
		Posit/ApproxDiv/PositApproxDiv.cpp \
		Posit/LogApproxMult/PositLAM.cpp \
		Posit/Cmp/PositComparator.cpp \
		Posit/PositAssign.cpp

libflopoco_la_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/FPExpressions

VHDLLexer.cpp:$(top_srcdir)/src/VHDLLexerRule.l
	${LEX} -o VHDLLexer.cpp $(top_srcdir)/src/VHDLLexerRule.l

VHDLLexer.hpp:VHDLLexer.cpp

Factories.cpp:$(top_srcdir)/src/Factories/MakeFactory.sh
	$(top_srcdir)/src/Factories/MakeFactory.sh $(top_srcdir)/src/Factories $(abs_builddir)/Factories.cpp

ExpressionParser.cpp: $(top_srcdir)/src/FPExpressions/ExpressionParser.y ExpressionScanner.cpp
	${YACC} -y -d --defines=$(top_srcdir)/src/FPExpressions/ExpressionParser.h --output=ExpressionParser.cpp --name-prefix=FlopocoExpression $(top_srcdir)/src/FPExpressions/ExpressionParser.y

ExpressionScanner.cpp:$(top_srcdir)/src/FPExpressions/ExpressionScanner.l
	${LEX} --prefix=FlopocoExpression --header-file=$(top_srcdir)/src/FPExpressions/FlopocoScanner.h -o ExpressionScanner.cpp  $(top_srcdir)/src/FPExpressions/ExpressionScanner.l

$(top_srcdir)/src/FPExpressions/FlopocoScanner.h:ExpressionScanner.cpp ExpressionParser.cpp

$(top_srcdir)/src/FPExpressions/ExpressionParser.h:ExpressionParser.cpp
