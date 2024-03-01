
# See output in <BUILD_DIR>/Testing/Temporary/LastTest.log

foreach(test_name
  atom_coord
  constraint_1
  constraint_2
  constraint_3
  disp_1a
  disp_1b
  disp_1c
  disp_2a
  disp_2b
  disp_3a
  disp_3b
  efield
  elec_1a
  elec_1b
  elec_1c
  elec_2a
  elec_2b
  elec_3a
  elec_3b
  elpot
  grad_1
  hess_1
  md_1
  md_2
  md_3
  opt_1
  pairwise_0
  pairwise_1
  pairwise_2
  pairwise_x
  pbc_1
  pbc_2
  pol_1a
  pol_1b
  pol_2a
  pol_2b
  pol_3a
  pol_3b
  print
  qm_1a
  qm_1b
  qm_2a
  qm_2b
  reduced
  sp_1
  symm_1
  symm_2full
  symm_2
  symm_2pw
  symm_2pw_printout
  total_1a
  total_2a
  total_3a
  total_4a
  total_4b
  total_4c
  total_4d
  total_5a
  total_5b
  total_5c
  total_5d
  total_6a
  total_6b
  total_6c
  total_6d
  xr_1
  xr_2
  xr_3
  # crambin/opt.in not handled (cannot load ff geometry)
  )

    add_test(
      NAME ${test_name}
      COMMAND efpmd ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.in
      )
    set_tests_properties(
      ${test_name}
      PROPERTIES
        PASS_REGULAR_EXPRESSION "COMPLETED SUCCESSFULLY"
        FAIL_REGULAR_EXPRESSION "DOES NOT MATCH"
      )
endforeach()
