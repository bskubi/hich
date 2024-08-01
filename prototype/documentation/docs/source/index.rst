.. Hich documentation main file

.. image:: media/hich_banner.png


Features:
   - **Comprehensive** .fastq/.bam/.pairs -> .hic/.mcool -> compartments, insulation, loops
   - **Merge** techreps -> bioreps -> conditions, matching coverage for fair feature calling.
   - **QC** data with a `read-level MultiQC report <https://multiqc.info/example-reports/hi-c/>`_ and `Hicrep <https://genome.cshlp.org/content/27/11/1939.short>`_ matrix-vs-matrix comparisons.
   - **Auto-setup** needed resources, such as genome references, indexes, and fragment files.

..
      Data from
            Paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8446342/
            4DN data page: https://data.4dnucleome.org/publications/dfc530f1-82c0-4ddc-8f95-6f40417f87a0/
      Basic workflow control (stub, humid, last step, resume)
      Set up the samples.csv spreadsheet
      Modify nextflow.config
      Full vignette with the comparison between different assay parameters


.. toctree::

      basicworkflowcontrol

.. note::

      This project is under active development.