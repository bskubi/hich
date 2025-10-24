package hich
import groovy.text.SimpleTemplateEngine

trait TaskPlan {
    Map<String, Command> bind_script = [:]
    Map<String, Object> bind_stub = [:]
    String script_template
    String stub_template

    String getScript() {
        def engine = new SimpleTemplateEngine()
        return engine.createTemplate(this.script_template).make(this.bind_script).toString()
    }

    String getStub() {
        def engine = new SimpleTemplateEngine()
        return engine.createTemplate(this.stub_template).make(this.bind_stub).toString()
    }
}