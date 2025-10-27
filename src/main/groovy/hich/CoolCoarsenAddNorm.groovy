package hich

class CoolCoarsenAddNorm extends TaskPlan {
    String mcool

    CoolCoarsenAddNorm ( 
        id, 
        cool,
        config_cool_coarsen_addnorm,
        cpus
    ) {
        this.script_template = '${command_cooler_zoomify}'
        this.stub_template = 'touch \'${mcool}\''
        this.mcool = "${id}.mcool"

        def bind = [
            id: id, 
            cool: cool.toString(),
            mcool: this.mcool,
            cpus: cpus
        ]

        this.bind_script.command_cooler_zoomify = new Command(config_cool_coarsen_addnorm.command_cooler_zoomify, bind)
        this.bind_stub = [mcool: this.mcool]
    }
}