#include <gtest/gtest.h>

// Include the header file that contains the compute_L2_norm function
#include "minimal_surface.h"

// Define a fixture class for the tests
class MinimalBenchmarkTest : public ::testing::Test {
protected:
  void SetUp() override {
    // Set up any necessary objects or variables for the tests
  }

  void TearDown() override {
    // Clean up any resources used by the tests
  }

  // Define any helper functions that you may need for the tests
};

// Test case for the compute_L2_norm function
TEST_F(MinimalBenchmarkTest, ComputeL2NormTest) {
  // Create an instance of the Minimal_Benchmark class
  Minimal_Benchmark minimalBenchmark;

  // Create a Vector<double> object for testing
  Vector<double> vec;

  // Call the compute_L2_norm function and get the result
  double result = minimalBenchmark.compute_L2_norm(vec);

  // Define the expected result
  double expected = 0.0; // Replace with the expected result

  // Check if the result matches the expected result
  EXPECT_DOUBLE_EQ(result, expected);
}

// Run all the tests
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}