add_custom_target(runVerificationTest 
  OUTPUT verification.out
  COMMAND echo $<TARGET_FILE:OpenSeesTcl> runVerificationSuite.tcl)

add_test(NAME VerificationTest COMMAND ${CMAKE_COMMAND} --build . --target runVerificationTest)

