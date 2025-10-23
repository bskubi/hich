import groovy.text.SimpleTemplateEngine

process FASTQ_ALIGN {
    container 'bskubi/hich_alignment:latest'
    debug true

    input:
    tuple val(id), 
          path(fastq, arity: '1..2'), 
          val(aligner), 
          path(aligner_index_dir), 
          val(aligner_index_prefix),
          val(config)

    output:
    tuple val(id), path(bam), emit: bam

    script:
    bam = "${id}.bam"
    fastq = fastq
        .withIndex()
        .collectEntries{f, i -> ["fastq${i+1}".toString(), f]}
    bind = [
        id: id, 
        *: fastq, 
        aligner: aligner,
        aligner_index_dir: aligner_index_dir, 
        aligner_index_prefix: aligner_index_prefix,
        'cpus': "${task.cpus}",
        bam: bam
    ]
    def engine = new SimpleTemplateEngine()
    Map align_command = config.align_command

    base_command = engine.createTemplate(align_command.base_command).make(bind).toString()
    flags = align_command.flags.findAll().collect {flag ->
        engine.createTemplate(flag).make(bind).toString()
    }
    opts = align_command.opts.findAll{ k, v -> v }.collect { k, v ->
        k_fmt = engine.createTemplate(k).make(bind).toString()
        v_fmt = engine.createTemplate(v).make(bind).toString()
        "${k_fmt} '${v_fmt}'"
    }
    args = align_command.args.findAll().collect { arg ->
        engine.createTemplate(arg).make(bind).toString()
    }

    align_command_script = [base_command, *flags, *opts, *args].join(" ")
    println(align_command_script)
    align_command_script
    

    stub:
    bam = "${id}.bam"
    """
    touch '${bam}'
    """
}