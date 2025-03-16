Install Hich and its dependencies
=================================

This guide shows how to install Hich's dependencies using `mamba <https://mamba.readthedocs.io/en/latest/index.html>`_.

1. `Nextflow <https://nextflow.io/docs/latest/index.html>`_
2. A `container management system <https://en.wikipedia.org/wiki/Containerization_(computing)>`_ like `Docker <https://docs.docker.com/>`_ or, for shared computing environments, `Singularity <https://docs.sylabs.io/guides/latest/user-guide/>`_/`Apptainer <https://apptainer.org/>`_.
3. `Squashfuse <https://github.com/vasi/squashfuse>`_, an up-to-date version of `OpenJDK <https://openjdk.org/>`_, `pip <https://pip.pypa.io/en/stable/>`_ and the `Hich CLI utilities <https://pypi.org/project/hich/>`_.

.. note::
    Nextflow, and therefore Hich, is primarily designed to be run on Linux or macOS, but can be installed on Windows via WSL. These instructions have been tested on Linux and may need adaptation to other operating systems. See `Nextflow's installation instructions <https://www.nextflow.io/docs/latest/install.html>`_ for details.

.. note::
    Assuming most users will run Hich on a shared computing environment, this documentation uses singularity. However, `Nextflow supports many container management systems <https://nextflow.io/docs/latest/container.html>`_ and it should be possible to use any of them Hich as long as it's compatible with your computing environment.

.. note::
    Hich uses the container management system to install further dependencies the first time you run it. It can be configured to permanently store these additional dependencies in a centralized location with the ``.cacheDir`` setting in the ``nextflow.config`` file. See `Nextflow's container documentation <https://nextflow.io/docs/latest/container.html>`_ for more details.

-----------------------------------------------------------------------------------------------------

Create the ``hich`` environment and install dependencies
........................................................


.. code-block:: c

    mamba create -qy -n hich bioconda::nextflow conda-forge::singularity conda-forge::squashfuse conda-forge::openjdk=21 pip && mamba activate hich && pip install hich

You can activate this environment with ``mamba activate hich``

.. note::
    ``mamba`` and ``conda`` can be used interchangeably, but ``mamba`` is typically faster.

Download the Hich pipeline
...........................

Create a directory where the ``hich`` pipeline code will be stored, go there, then clone it from github.

.. code-block:: c

    mkdir ~/hich && cd ~/hich
    git clone https://github.com/bskubi/hich && cd hich
