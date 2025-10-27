package hich
import groovy.text.SimpleTemplateEngine

class TaskPlan {
    Map<String, Command> bind_script = [:]
    Map<String, Object> bind_stub = [:]
    String script_template
    String stub_template

    String getScript() {
        SimpleTemplateEngine engine = new SimpleTemplateEngine()
        return engine.createTemplate(this.script_template).make(this.bind_script).toString()
    }

    String getStub() {
        SimpleTemplateEngine engine = new SimpleTemplateEngine()
        return engine.createTemplate(this.stub_template).make(this.bind_stub).toString()
    }

    String toString() {
        Map parts = [
            script: getScript(),
            stub: getStub()
        ]
        return parts.toString()
    }
}