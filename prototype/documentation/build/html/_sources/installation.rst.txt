Install Hich
============

`Install with Mamba <https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html>`_
-----------------------------------------------------------------------------------------------------

Install

.. code-block:: c

    mamba create -qy -c bioconda -n hich nextflow singularity

Activate environment (each time you log in, prior to running Hich)

.. code-block::

    mamba activate hich

`Install with Conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_
-----------------------------------------------------------------------------------------------------

Install

.. code-block:: c

    conda create -qy -c bioconda -n hich nextflow singularity

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
    | ``-c bioconda`` lets mamba/conda look for packages in the bioconda repository, which is where the nextflow package is located.
    | ``-n hich`` names the environment "hich". You can choose a different name.
    | ``nextflow singularity`` The two packages that must be installed to the ``hich`` environment.

[mamba|conda] activate
.......................
Each time you log in, Mamba/Conda load a base or default environment. Each time you log in, activate this environment with ``mamba activate hich`` or ``conda activate hich``. Then you can run Hich.

The outputs will be accessible to you whether or not the environment is activated.

New to Nextflow and Singularity?
--------------------------------

`Nextflow <https://nextflow.io/>`_ is a powerful workflow platform for reproducible science, similar to `Snakemake <https://snakemake.github.io/>`_ and widely used in bioinformatics. `Singularity <https://docs.sylabs.io/guides/3.5/user-guide/introduction.html>`_ makes software extremely portable. The ``nextflow run bskubi/hich`` command downloads ``hich`` from its `GitHub repo <https://github.com/bskubi/hich>`_, installs it to the Nextflow home directory, and runs it using the :doc:`sample file <samplefile>` in your current working directory with default settings. See Nextflow's section on `pipeline sharing <https://www.nextflow.io/docs/latest/sharing.html>`_ for more details.

For Hich's individual processes, from data processing to QC to analysis, Hich relies on a number of existing and novel tools, packages and algorithms. These have been containerized for Hich so that Singularity can install them behind the scenes, without any need for manual setup or configuration by the user.