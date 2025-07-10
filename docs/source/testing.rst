
=======
Testing
=======



A set of tests are present in the ``tests-integration`` folder. These cases are very simple and have no physical accuracy, are just meant to test the integrated working of the scripts functionalities.
To run these execute the following commands from a bash shell inside the tests folder called ``FLUNED-Repository/tests-integration``.

1. Run the test suite.

.. code-block:: bash
   :linenos:

    ./run_tests.sh

2. Run the post processing script for the tests simulations:

.. code-block:: bash
   :linenos:

    ./run_post.sh

3. Run a script that executes the comparison between the test results and the expected results. This comparison will be written in the ``test_results`` file
   
.. code-block:: bash
   :linenos:

    ./eval_results.sh

4. Clean the created files

.. code-block:: bash
   :linenos:

    ./clean_tests.sh

