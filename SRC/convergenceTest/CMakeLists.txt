add_library(convergenceTest "")

target_sources(convergenceTest
    PRIVATE
        ConvergenceTest.cpp
	CTestNormUnbalance.cpp
	CTestNormDispIncr.cpp
	CTestEnergyIncr.cpp
	CTestRelativeNormUnbalance.cpp
	CTestRelativeNormDispIncr.cpp
	CTestRelativeEnergyIncr.cpp
	CTestFixedNumIter.cpp
	CTestRelativeTotalNormDispIncr.cpp
	NormDispAndUnbalance.cpp
	NormDispOrUnbalance.cpp
	CTestPFEM.cpp
    PUBLIC
        ConvergenceTest.h
	CTestNormUnbalance.h
	CTestNormDispIncr.h
	CTestEnergyIncr.h
	CTestRelativeNormUnbalance.h
	CTestRelativeNormDispIncr.h
	CTestRelativeEnergyIncr.h
	CTestFixedNumIter.h
	CTestRelativeTotalNormDispIncr.h
	NormDispAndUnbalance.h
	NormDispOrUnbalance.h
	CTestPFEM.h
)

target_include_directories(convergenceTest PUBLIC $(CMAKE_CURRENT_LIST_DIR))
