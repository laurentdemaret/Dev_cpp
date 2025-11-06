
add_executable (TestCurveFitting curve_fitting.cc )
target_link_libraries(TestCurveFitting Eigen3::Eigen ${CERES_LIBRARIES} ) 
