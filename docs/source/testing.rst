=======
Testing
=======

A set of tests are present in the ``tests-integration`` folder. 
These cases are very simple and have no physical accuracy; they are just meant to test the integrated working of the scripts' functionalities. 
To run them, execute the following commands from a Bash shell inside the ``FLUNED-Repository/tests-integration`` folder.

1. Run the test suite.

   .. code-block:: bash

      ./run_tests.sh

2. Run the post-processing script for the test simulations:

   .. code-block:: bash

      ./run_post.sh

3. Compare the test results with the expected results (output written to ``test_results``):

   .. code-block:: bash

      ./eval_results.sh

4. Clean the generated files.

   .. code-block:: bash

      ./clean_tests.sh
