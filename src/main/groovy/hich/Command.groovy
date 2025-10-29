package hich
import groovy.text.SimpleTemplateEngine
import java.util.Collection

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

        
        this.opts.findAll{ k, v -> v }
        .each { k, v ->
            /** Allow multi-options (i.e. -k val1 -k val2)

            Convert non-collection values to a list. Iterate through list
            of values and add repeat option for each element.
            */
            String k_fmt = engine.createTemplate(k).make(bind).toString()
            v = Collection.isInstance(v) ? v : [v]
            v.each { v_i ->
                String v_fmt = engine.createTemplate(v_i).make(bind).toString()
                command += ["${k_fmt} '${v_fmt}'"]
            }
        }

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