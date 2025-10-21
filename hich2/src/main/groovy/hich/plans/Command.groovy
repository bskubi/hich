package hich.plans
import hich.format.Format

class Command {
    String base_command = null
    Map opts = [:]
    Map args = [:]
    List prepared_opts = null
    List prepared_args = null
    String prepared_command = null

    Command setBaseCommand(String command) {
        base_command = command
        return this
    }

    Command updateOpts(Map new_opts) {
        if (new_opts) {
            opts += new_opts
        }
        
        return this
    }

    Command updateArgs(Map new_args) {
        if (new_args) {
            args += new_args
        }
        return this
    }

    Command formatBaseCommand(Map bind) {
        base_command = Format.format(base_command, bind)
        return this
    }

    Command formatOpt(String opt, Map bind) {
        opts = Format.updateFormat(opts, opt, bind)
        return this
    }

    Command formatArg(String arg, Map bind) {
        args = Format.updateFormat(args, arg, bind)
        return this
    }

    Command evalOpt(String opt, Map bind) {
        if (opts.containsKey(opt)) {
            bind.each {symbol, object ->
                opts[opt] = Eval.me(symbol, object, opts[opt])
            }
        }
        return this
    }

    Command evalArg(String arg, Map bind) {
        if (args.containsKey(arg)) {
            bind.each {symbol, object ->
                args[arg] = Eval.me(symbol, object, args[arg])
            }
        }
        return this
    }

    Command removeUnusedOpts() {
        opts = opts.findAll{it.value != null && it.value != false}
        return this
    }

    Command removeUnusedArgs() {
        args = args.findAll{it.value != null}
        return this
    }

    Command prepareOpts() {
        prepared_opts = opts.collect{k, v -> v == true ? k.toString() : "${k} '${v}'".toString()}
        return this
    }

    Command prepareArgs() {
        prepared_args = args.collect{k, v -> "'${v}'".toString()}
        return this
    }

    Command prepareCommand() {
        prepared_command = ([base_command] + prepared_opts + prepared_args).join(" ")
        return this
    }

    Command buildCommand() {
        return (
            removeUnusedOpts()
            .prepareOpts()
            .removeUnusedArgs()
            .prepareArgs()
            .prepareCommand()
        )
    }
}