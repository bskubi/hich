package hich
import groovy.text.SimpleTemplateEngine

class Command {
    String base_command
    Map opts
    List flags
    List args
    Map bind

    Command(Map config, Map bind) {
        this.base_command = config.base_command
        this.opts = config.opts
        this.flags = config.flags
        this.bind = bind
    }

    String toString() {
        def engine = new SimpleTemplateEngine()
        def command = []

        command.add(engine.createTemplate(base_command).make(bind).toString())

        command += (
            opts.findAll{ k, v -> v }
            .collect { k, v ->
                def k_fmt = engine.createTemplate(k).make(bind).toString()
                def v_fmt = engine.createTemplate(v).make(bind).toString()
                "${k_fmt} '${v_fmt}'"
            }
        )

        command += (
            flags.findAll()
            .collect { flag ->
                engine.createTemplate(flag).make(bind).toString()
            }
        )

        command += (
            args.findAll()
            .collect { arg ->
                engine.createTemplate(arg).make(bind).toString()
            }
        )

        return command.join(" ")
    }
}