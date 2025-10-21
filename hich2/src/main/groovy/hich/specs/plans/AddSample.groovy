package hich.specs.plans
import hich.specs.plans.keys.PlanKey

class AddSample {
    private static final Map spec = [
        (PlanKey.TRANSFORMERS): [
            joinPermissive: [
                [
                    output_key: "id", 
                    input_keys: ["condition", "biorep", "techrep"], 
                    separator: "_"
                ]
            ]
        ],
        (PlanKey.POST_VALIDATORS): [
            hasNonWhitespaceChar: [
                [key: "id"],
                [key: "assembly"]
            ]
        ]
    ]
}