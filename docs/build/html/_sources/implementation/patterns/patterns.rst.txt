Common Patterns
---------------------------

One -> One
....................................................

**Use case: processes that take individual samples as input and return new or updated attributes.**

**Filter** for the required samples, **extract** the sample attributes needed (usually including ``id``), call the **process**, call the process, **label** the outputs, and **pack** the outputs as new attributes of the same by using ``id`` to link the new outputs to the correct sample.

Here is an example from ``parse.nf``. 

.. code-block:: c

    samples
        | filter{!skip("parse") && it.datatype in ["fastq", "sambam"]}
        | map{tuple(it.id, it.sambam, it.chromsizes, it.assembly, it.pairtoolsParse2Params, it.reshapeParams, it.subMap("minMapq"))}
        | PairtoolsParse2
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    pack(samples, result) | set{samples}

Many -> One
.....................................................

``groupRowsToColumnFormat`` groups the samples by common values of ``biorep``, ``condition`` and ``aggregationPlanName`` to yield a channel with one column-format item per group. Then ``coalesce`` replaces any constant vectors with a constant value (i.e. ``[val: [1, 1, 1]]`` becomes ``[val: 1]``). 

Example:

.. code-block:: c

    // Filter samples
    levelSamples
    | branch {
            yesMerge: it.includeInMerge && it.mergeTechrepsToBioreps
            noMerge: true
        }
    | set{samples}

    // Merge the pairs files.
    groupRowsToColumnFormat(samples.yesMerge, ["biorep", "condition", "aggregationPlanName"], ["dropNull": true])
        | map{coalesce(it)}
        | map{tuple(makeID(it, columns = true), it.latestPairs)}
        | mergeTechrepsToBioreps
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{mergedTechrepAttributes}

    // Group the merged result by the mergeGroupIdentifiers, then coalesce common values
    // to a single value, dropping any null or heterogeneous values. The other
    // common values are kept as inherited merge attributes. Then add an ID
    // for the merge.
    groupRowsToColumnFormat(samples.yesMerge, ["biorep", "condition", "aggregationPlanName"], ["dropNull": true])
        | map{coalesce(it, '_drop')}
        | map{it += [id:makeID(it)]}
        | set{inheritedMergeAttributes}