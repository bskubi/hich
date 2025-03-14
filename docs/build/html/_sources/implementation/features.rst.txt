Sample Selection Strategy
=========================

strategyKeys

.. code-block:: c

  hicrep:
    hicrep_default:
      sampleSelectionStrategy: []
      resolutions: [10000]
      chroms: ["chr1", "chr2", "chr3"]
      exclude: ["chr3"]
      chromFilter: []
      h: [1]
      dBPMax: [5000000]
      bDownSample: [true]

createCompositeStrategy

A composite strategy is a hashmap in which keys are sample attributes and values are lists of permitted sample attribute values. It is created by combining one or more individual strategies specified in params.sampleSelectionStrategies.

filterSamplesByStrategy

After a composite strategy is built, filter for samples for which all sample attributes are present and are in the list of permitted values specified by the composite strategy.

pairSamplesByStrategy

Get all pairs of samples having matching values of strategy.same.

groupSamplesByStrategy

Get all samples having matching values of strategy.same. 

.. toctree::
    :hidden:

    implementation/index