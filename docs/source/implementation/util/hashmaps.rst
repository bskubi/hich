``util/hashmaps.nf``
------------------------------------

.. code-block:: c

    def asHashSet(val) {
        // Convert a non-HashSet val into a HashSet
        if (val instanceof ArrayList) return new HashSet(val)
        if (val instanceof HashSet) return val
        return new HashSet([val])
    }