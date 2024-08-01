How to launch Hich stub, humid and full runs
==========================================================================

.. note::
    This vignette shows how to run Hich using any of four dependency management solutions:
    Singularity, Docker, Conda or Mamba. Choose the version that fits your setup.


1. Install `Nextflow <https://www.nextflow.io/docs/latest/install.html>`_ and one of the following: `Singularity <https://docs.sylabs.io/guides/3.0/user-guide/installation.html>`_, `Docker <https://docs.docker.com/engine/install/>`_, `Conda <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_, or `Mamba <https://www.google.com/search?q=mamba+install&oq=mamba+install&gs_lcrp=EgZjaHJvbWUyFAgAEEUYORhDGIMBGLEDGIAEGIoFMgwIARAAGEMYgAQYigUyDAgCEAAYQxiABBiKBTIMCAMQABhDGIAEGIoFMgwIBBAAGEMYgAQYigUyBggFEEUYPDIGCAYQRRg8MgYIBxBFGDzSAQc5NjZqMGo3qAIAsAIA&sourceid=chrome&ie=UTF-8>`_

2. Clone Hich into your home directory.

``git clone https://github.com/bskubi/hich.git ~ && cd ~/hich``

3. The vignette data files and experiment description are located in
``hich/vignettes``. Start with a stub run, which quickly verifies the workflow
is able to run without doing any data processing.

| ``nextflow run hich.nf -stub-run -singularity.enabled true``
| ``nextflow run hich.nf -stub-run -docker.enabled true``
| ``nextflow run hich.nf -stub-run -conda.enabled true``
| ``nextflow run hich.nf -stub-run -conda.enabled true -conda.useMamba true``

4. Next, do a "humid" run, which auto-downsamples the input data to the first 100,00 reads per file to quickly verify
Hich can process the data successfully.

| ``nextflow run hich.nf --humid -singularity.enabled true``
| ``nextflow run hich.nf --humid -docker.enabled true``
| ``nextflow run hich.nf --humid -conda.enabled true``
| ``nextflow run hich.nf --humid -conda.enabled true -conda.useMamba true``

5. Finally, process the entire dataset.

| ``nextflow run hich.nf -singularity.enabled true``
| ``nextflow run hich.nf -docker.enabled true``
| ``nextflow run hich.nf -conda.enabled true``
| ``nextflow run hich.nf -conda.enabled true -conda.useMamba true``