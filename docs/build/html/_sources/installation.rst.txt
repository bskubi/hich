Install Hich
============

`Install with Mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
-----------------------------------------------------------------------------------------------------

Install

.. code-block:: c

    mamba create -qy -n hich bioconda::nextflow conda-forge::squashfuse conda-forge::apptainer

Activate environment (each time you log in, prior to running Hich)

.. code-block::

    mamba activate hich

`Install with Conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
-----------------------------------------------------------------------------------------------------

Install

.. code-block:: c

    conda create -qy -n hich bioconda::nextflow conda-forge::squashfuse conda-forge::apptainer

Activate environment (each time you log in, prior to running Hich)

.. code-block::

    conda activate hich

Run Hich
--------

First, set up your :doc:`sample file <samplefile>` in the current working directory. Then run Hich:

.. code-block:: c
    
    nextflow run bskubi/hich

New to Mamba/Conda?
-------------------

Mamba and Conda make it easy to install new packages. Mamba is a faster drop-in replacement for Conda. Here is an explanation of what these commands do:

[mamba|conda] create
.....................

The ``create`` command creates a new environment especially for Hich with its dependencies.
    | ``-qy`` quiets installation messages automatically says "yes" to the "do you want to install" prompt.
    | ``-n hich`` names the environment "hich". You can choose a different name.
    | ``bioconda::nextflow conda-forge::squashfuse conda-forge::apptainer`` The three packages that must be installed to the ``hich`` environment along with the channels from which they should be accessed. Note that it is essential to use these specific channels.

[mamba|conda] activate
.......................
Each time you log in, Mamba/Conda load a base or default environment. Each time you log in, activate this environment with ``mamba activate hich`` or ``conda activate hich``. Then you can run Hich.

The outputs will be accessible to you whether or not the environment is activated.

New to Nextflow, Squashfuse or Apptainer?
-------------------------------------------

`Nextflow <https://nextflow.io/>`_ is a powerful workflow platform for reproducible science, similar to `Snakemake <https://snakemake.github.io/>`_ and widely used in bioinformatics. `Apptainer <https://apptainer.org/docs/admin/latest/index.html>`_ makes software extremely portable and works in a variety of HPC environments. `Squashfuse <https://github.com/vasi/squashfuse>`_ is a critical to Apptainer's performance. The ``nextflow run bskubi/hich`` command downloads ``hich`` from its `GitHub repo <https://github.com/bskubi/hich>`_, installs it to the Nextflow home directory, and runs it using the :doc:`sample file <samplefile>` in your current working directory with default settings. See Nextflow's section on `pipeline sharing <https://www.nextflow.io/docs/latest/sharing.html>`_ for more details.

For Hich's individual processes, from data processing to QC to analysis, Hich relies on a number of existing and novel tools, packages and algorithms. These have been containerized for Hich so that Apptainer can install them behind the scenes, without any need for manual setup or configuration by the user.