def validateMemory(memory, buffer, minimum) {
    // Use 2G less than the total memory allocated for the job
    // to a minimum of 2G
    memory ? Math.max(memory.toGiga() - buffer, minimum) : minimum
}
