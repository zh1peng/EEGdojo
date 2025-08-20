# EEGdojo Toolbox Test Report

## Overall Summary

All four test scripts failed to run to completion. The failures point to two main categories of issues:

1.  **Missing External Dependencies**: The tests that use the `remove_powerline` function are failing because the `pop_cleanline` function, which is part of the EEGLAB toolbox, could not be found.
2.  **MATLAB Path and Private Functions**: The tests that use the `remove_bad_channels` function are failing because the helper function `cleanraw_rejchan`, which is located in the `+EEGdojo/private` folder, could not be found. This indicates an issue with how MATLAB is resolving the path for package-private functions when the tests are executed.

Below is a detailed analysis of each test.

---

## Test-by-Test Analysis

### 1. Test Script: `run_all_tests_simple.m`

*   **Result**: FAILED
*   **Error**: `Undefined function 'pop_cleanline' for input arguments of type 'struct'.`
*   **Analysis**: The script executed several tests successfully but failed when it reached the `EEGdojo.preprocess.remove_powerline` function. This function internally calls `pop_cleanline`, which is a standard function from the EEGLAB toolbox.
*   **Suggested Fix**: This is an environment issue. Please ensure that your EEGLAB installation is correctly added to your MATLAB path before running the tests. You can do this with the `addpath` command, for example: `addpath('C:\path\to\eeglab')`.

### 2. Test Script: `test_pipeline.m`

*   **Result**: FAILED
*   **Error**: `Undefined function 'pop_cleanline' for input arguments of type 'struct'.`
*   **Analysis**: This test failed for the same reason as the previous one. The pipeline was executing correctly until it reached the `remove_powerline` step, at which point it failed due to the missing `pop_cleanline` function.
*   **Suggested Fix**: The same fix applies here: ensure EEGLAB is on the MATLAB path.

### 3. Test Script: `test_bad_channel_detection.m`

*   **Result**: FAILED
*   **Error**: `Undefined function 'cleanraw_rejchan' for input arguments of type 'struct'.`
*   **Analysis**: This test failed inside the `EEGdojo.preprocess.remove_bad_channels` function. This function attempts to call `cleanraw_rejchan`, which is located in the `+EEGdojo/private` directory. Although the test script adds the project's root directory to the path, MATLAB is still unable to find this private function. This can happen if the context from which the test is run does not correctly inform MATLAB about the package structure.
*   **Suggested Fix**: The most robust solution is to have a single test runner script in the project root that sets up the path correctly for all tests. While I have provided scripts that attempt to do this, MATLAB's handling of package-private functions can be sensitive. You may need to ensure that you are running the tests from the project root directory, or that the `+EEGdojo` package is fully recognized by MATLAB by having its parent directory on the path at startup.

### 4. Test Script: `test_pipeline_with_bad_channels.m`

*   **Result**: FAILED
*   **Error**: `Undefined function 'cleanraw_rejchan' for input arguments of type 'struct'.`
*   **Analysis**: I initially fixed a bug in this test related to how parameters were being passed to the `Params` class. After the fix, the test now fails for the same reason as `test_bad_channel_detection.m`. The `remove_bad_channels` function, when called from within the pipeline, cannot locate its private helper function `cleanraw_rejchan`.
*   **Suggested Fix**: The same fix as for the previous test applies. The core issue is the accessibility of private package functions in your test environment.
