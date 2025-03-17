Install Hich and its dependencies
=================================

This guide shows how to install Hich's dependencies using `mamba <https://mamba.readthedocs.io/en/latest/index.html>`_. 

-----------------------------------------------------------------------------------------------------

Create the ``hich`` environment and install dependencies
........................................................


.. code-block:: c

    mamba create -qy -n hich bioconda::nextflow conda-forge::singularity conda-forge::squashfuse conda-forge::openjdk=21 pip && mamba activate hich && pip install hich

The dependencies that will be installed by this command are:

+  `Nextflow <https://nextflow.io/docs/latest/index.html>`_
+  `Singularity <https://docs.sylabs.io/guides/latest/user-guide/>`_, which is a `container management system <https://en.wikipedia.org/wiki/Containerization_(computing)>`_ similar to `Docker <https://docs.docker.com/>`_ but more commonly compatible with shared computing environments
+  `Squashfuse <https://github.com/vasi/squashfuse>`_, an up-to-date version of `OpenJDK <https://openjdk.org/>`_, `pip <https://pip.pypa.io/en/stable/>`_ and the `Hich CLI utilities <https://pypi.org/project/hich/>`_


This will activate the ``hich`` environment. You will need to activate it prior to running the pipeline via ``mamba activate hich``.

.. note::
    ``mamba`` and ``conda`` can be used interchangeably, but ``mamba`` is typically faster.

Download the Hich pipeline
...........................

Create a directory where the ``hich`` pipeline code will be stored, go there, then clone it from github.

.. code-block:: c

    mkdir ~/hich && cd ~/hich
    git clone https://github.com/bskubi/hich && cd hich

.. note::
    Nextflow, and therefore Hich, is primarily designed to be run on Linux or macOS, but can be installed on Windows via WSL. These instructions have been tested on Linux and may need adaptation to other operating systems. See `Nextflow's installation instructions <https://www.nextflow.io/docs/latest/install.html>`_ for details.

.. note::
    Assuming most users will run Hich on a shared computing environment, this documentation uses singularity. However, `Nextflow supports many container management systems <https://nextflow.io/docs/latest/container.html>`_ and it should be possible to use any of them Hich as long as it's compatible with your computing environment.

.. note::
    Hich uses the container management system (i.e. Singularity) to install further dependencies the first time you run it. It can be configured to permanently store these additional dependencies in a centralized location with the ``.cacheDir`` setting in the ``nextflow.config`` file. See `Nextflow's container documentation <https://nextflow.io/docs/latest/container.html>`_ for more details.