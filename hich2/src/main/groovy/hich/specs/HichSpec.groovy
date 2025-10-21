package hich.specs
import hich.specs.plans.Align
import hich.specs.plans.AddSample
import hich.specs.plans.keys.PlanKey



class HichSpec {
    enum SpecKey {
        ADD_SAMPLE,
        ALIGN
    }

    static {
        def all_defined_keys = PlanKey.values().collect { it.toString() }.toSet()

        SPECS.each { name, spec ->
            def defined_keys = spec.keySet().collect { it.toString() }.toSet()
            def undefined_keys = defined_keys - all_defined_keys
            if (!undefined_keys.isEmpty()) {
                throw new Exception("FATAL: HichSpec for spec '${name}' has undefined keys: ${undefined_keys}")
            }
        }
    }

    private static final Map SPECS = [
        (SpecKey.ADD_SAMPLE): AddSample.spec,
        (SpecKey.ALIGN): Align.spec
    ]
}