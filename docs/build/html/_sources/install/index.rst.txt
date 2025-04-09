Install Hich and its dependencies
=================================

This guide shows how to install Hich's dependencies using `mamba <https://mamba.readthedocs.io/en/latest/index.html>`_. Hich has two components: a Nextflow pipeline and a Python CLI app. The pipeline uses the CLI app for some operations, and the CLI app also has useful operations independent of the pipeline.

Create the ``hich`` environment and install dependencies
........................................................

.. code-block:: c

    mamba env create -f pipeline_env.yml && mamba activate hich


This will activate the ``hich`` environment. You will need to activate again prior to running the pipeline with ``mamba activate hich``.

.. note::
    ``mamba`` and ``conda`` can be used interchangeably, but ``mamba`` is typically faster.

Download the Hich pipeline
...........................

.. code-block:: c

    git clone "https://github.com/bskubi/hich" && cd hich

.. note::
    Assuming most users will run Hich on a shared computing environment, this documentation assumes you are using Singularity as your containerization solution when running Nextflow. However, `Nextflow supports many container management systems <https://nextflow.io/docs/latest/container.html>`_ and it should be possible to use any of them Hich as long as it's compatible with your computing environment. Nextflow can be configured to cache containers with the ``.cacheDir`` setting in the ``nextflow.config`` file. See `Nextflow's container documentation <https://nextflow.io/docs/latest/container.html>`_ for more details.

.. note::
    Nextflow, and therefore Hich, is primarily designed to be run on Linux or macOS, but can be installed on Windows via WSL. These instructions have been tested on Linux and may need adaptation to other operating systems. See `Nextflow's installation instructions <https://www.nextflow.io/docs/latest/install.html>`_ for details.