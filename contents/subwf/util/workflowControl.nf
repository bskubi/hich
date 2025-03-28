/*
    emptyOnLastStep is used to control behavior when the workflow reaches the step specified by the --lastStep
    param, if specified. It should get passed the main samples channel and returns an empty channel if it's the 
    last step (halting execution) or the original samples channel otherwise.

    If --viewLastStep is specified it will display the contents of the samples channel after the last step finishes.
    If a space-separated string is specified, it will turn that into a list and submap just those sample attributes for viewing.
*/
def emptyOnLastStep(step, samples) {
    def isExplicitLastStep = (params.containsKey("lastStep") && params.get("lastStep") == step)
    def isLastStep = (step == "End") || isExplicitLastStep
    def hasViewLastStep = params.containsKey("viewLastStep") && params.get("viewLastStep")
    if (isLastStep && hasViewLastStep) {
        samples
            | map {
                sample ->
                
                params.viewLastStep instanceof Boolean ? sample : sample.subMap(params.viewLastStep.split())}
            | view
    }
    return isExplicitLastStep ? channel.empty() : samples
}

def skip(step) {
    /*
        Users may want to skip some steps, such as QC or forming a particular kind of contact matrix,
        or run only certain steps. This uses both params to define a list of steps to be skipped
        (the intersection of skip and runOnly's complement).
    */
    def excluded = params.containsKey("runOnly") && !params.runOnly.split().contains(step)
    def skipped = params.containsKey("skip") && params.skip.split().contains(step)
    return excluded || skipped
}
