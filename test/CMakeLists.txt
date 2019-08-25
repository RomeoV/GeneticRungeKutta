find_package(GTest REQUIRED)

add_executable(simpson_rule simpson_rule.cpp)
target_link_libraries(simpson_rule 
  GeneticScheme
  GTest::GTest
  GTest::Main
  fmt
)
set_target_properties(simpson_rule
  PROPERTIES
  CXX_STANDARD 17
)
gtest_discover_tests(simpson_rule)
