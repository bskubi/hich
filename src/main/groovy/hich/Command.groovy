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
        this.args = config.args
        this.bind = bind
    }

    String toString() {
        SimpleTemplateEngine engine = new SimpleTemplateEngine()
        List<String> command = []

        command += [engine.createTemplate(this.base_command).make(this.bind).toString()]

        command += (
            this.opts.findAll{ k, v -> v }
            .collect { k, v ->
                String k_fmt = engine.createTemplate(k).make(bind).toString()
                String v_fmt = engine.createTemplate(v).make(bind).toString()
                "${k_fmt} '${v_fmt}'"
            }
        )

        command += (
            this.flags.findAll()
            .collect { flag ->
                engine.createTemplate(flag).make(bind).toString()
            }
        )

        command += (
            this.args.findAll()
            .collect { arg ->
                engine.createTemplate(arg).make(bind).toString()
            }
        )

        return command.join(" ")
    }
}