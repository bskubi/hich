Install Nextflow and Apptainer
==============================

Hich depends on Nextflow and a containerization solution such as Apptainer, which are easiest to install on a macOS or Linux-based system that has conda or mamba preinstalled. Nextflow can also be installed on Windows using WSL. For a conda-based installation, replace `mamba` with `conda` in the installation steps below. It will install Nextflow, the needed JDK version, apptainer, and the squashfuse dependency needed to efficiently start container instances.  

`Install with Mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
-----------------------------------------------------------------------------------------------------

Install

.. code-block:: c

    mamba create -qy -n hich bioconda::nextflow conda-forge::squashfuse conda-forge::apptainer conda-forge::openjdk=21

Activate environment (each time you log in, prior to running Hich)

.. code-block::

    mamba activate hich

Download Hich and vignettes
...........................

.. code-block:: c

    git clone https://github.com/bskubi/hich && cd hich

Run end-to-end test vignettes
.............................

If running on a multi-user environment (i.e. SLURM), allocate resources and create a batch execution script if desired. The command to launch the automated test scripts depends on whether you wish to launch a resource-constrained test (localPC) or higher resource test (localHPC). 

.. code-block:: c

    bash tests/test_localPC.sh

.. code-block:: c

    bash tests/test_localHPC.sh


New to Nextflow, Squashfuse or Apptainer?
-----------------------------------------

`Nextflow <https://nextflow.io/>`_ is a powerful workflow platform for reproducible science, similar to `Snakemake <https://snakemake.github.io/>`_ and widely used in bioinformatics. `Apptainer <https://apptainer.org/docs/admin/latest/index.html>`_ makes software extremely portable and works in a variety of HPC environments. `Squashfuse <https://github.com/vasi/squashfuse>`_ is a critical to Apptainer's performance. The ``nextflow run bskubi/hich`` command downloads ``hich`` from its `GitHub repo <https://github.com/bskubi/hich>`_, installs it to the Nextflow home directory, and runs it using the :doc:`sample file <samplefile>` in your current working directory with default settings. See Nextflow's section on `pipeline sharing <https://www.nextflow.io/docs/latest/sharing.html>`_ for more details.

For Hich's individual processes, from data processing to QC to analysis, Hich relies on a number of existing and novel tools, packages and algorithms, known as the "Hich toolkit." This toolkit has been installed to a frozen container image from which they can be run by a container management system (Apptainer, Singularity, Podman, Docker, etc). This image is prebuilt and runs with no installation anywhere the container management system exists, eliminating the need for cumbersome dependency installation found in other processing workflows.
