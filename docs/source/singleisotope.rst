Single Isotope workflow
-----------------------

For a single-isotope computation, the simulation workflow requires a completed CFD simulation.
The path to the simulation is provided through the fluned pre-processor input.
The pre-processor will analyize the CFD simulation and create inside it a FLUNED case simulation folder, according to the various parameters included in the input.
The various parameters are described in more detail in: :doc:`keywords`.

In addition, the test simulations in the  `FLUNED-Repository/tests-integration` folder constitute a useful resource to understand and test the wrokflow for single-isotope calculations.

From the command line, in the folder where the Ansys‑FLUENT or OpenFOAM files are located, the following workflow for single-isotope simulations is recommended.

1. Generate a FLUNED input template:

.. code-block:: bash

    fluned -t

2. Modify the input_template file with the appropriate parameters. The file will contain the default values together with the available options.

3. Launch the FLUNED pre‑processor:

.. code-block:: bash

    fluned -i input-file

This command creates a folder named as the FLUNED case reported in the input file.
The generated simulation files follow the OpenFOAM structure, so additional tweaking of
simulation parameters is possible if you have previous knowledge of OpenFOAM.
Thanks to the modularity of OpenFOAM classes, FLUNED simulations can be run in parallel
using the standard OpenFOAM commands.

4. In the FLUNED case folder, launch the solver:

.. code-block:: bash

    FLUNED-solver

5. Once the simulation is finished, run the post‑processor in the simulation folder:

.. code-block:: bash

    fluned-post

This writes a text summary to ``RESULTS/SUMMARY.csv``.
For steady-state cases it also writes a machine-readable summary to
``RESULTS/fluned_summary.xml``.
