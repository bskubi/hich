class HichCommandBuilder {
    String base_command = null
    Map opts = [:]
    Map args = [:]
    List prepared_opts = null
    List prepared_args = null
    String prepared_command = null

    HichCommandBuilder setBaseCommand(String command) {
        base_command = command
        return this
    }

    HichCommandBuilder updateOpts(Map new_opts) {
        if (new_opts) {
            opts += new_opts
        }
        
        return this
    }

    HichCommandBuilder updateArgs(Map new_args) {
        if (new_args) {
            args += new_args
        }
        return this
    }

    HichCommandBuilder formatBaseCommand(Map bind) {
        base_command = HichUtil.format(base_command, bind)
        return this
    }

    HichCommandBuilder formatOpt(String opt, Map bind) {
        opts = HichUtil.formatIfPresent(opts, opt, bind)
        return this
    }

    HichCommandBuilder formatArg(String arg, Map bind) {
        args = HichUtil.formatIfPresent(args, arg, bind)
        return this
    }

    HichCommandBuilder evalOpt(String opt, Map bind) {
        if (opts.containsKey(opt)) {
            bind.each {symbol, object ->
                opts[opt] = Eval.me(symbol, object, opts[opt])
            }
        }
        return this
    }

    HichCommandBuilder evalArg(String arg, Map bind) {
        if (args.containsKey(arg)) {
            bind.each {symbol, object ->
                args[arg] = Eval.me(symbol, object, args[arg])
            }
        }
        return this
    }

    HichCommandBuilder removeUnusedOpts() {
        opts = opts.findAll{it.value != null && it.value != false}
        return this
    }

    HichCommandBuilder removeUnusedArgs() {
        args = args.findAll{it.value != null}
        return this
    }

    HichCommandBuilder prepareOpts() {
        prepared_opts = opts.collect{k, v -> v == true ? k.toString() : "${k} '${v}'".toString()}
        return this
    }

    HichCommandBuilder prepareArgs() {
        prepared_args = args.collect{k, v -> "'${v}'".toString()}
        return this
    }

    HichCommandBuilder prepareCommand() {
        prepared_command = ([base_command] + prepared_opts + prepared_args).join(" ")
        return this
    }

    HichCommandBuilder buildCommand() {
        return (
            removeUnusedOpts()
            .prepareOpts()
            .removeUnusedArgs()
            .prepareArgs()
            .prepareCommand()
        )
    }
}