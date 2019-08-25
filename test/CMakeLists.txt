find_package(GTest REQUIRED)

add_executable(simpson_rule simpson_rule.cpp)
target_link_libraries(simpson_rule 
  GeneticScheme
  GTest::GTest
  GTest::Main
  fmt
)
gtest_discover_tests(simpson_rule)
