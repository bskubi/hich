package hich.plans
import java.nio.file.Path

class Plan {
    String command = null
    String stub = null
    Map<String, Object> input = null 
    Map<String, Object> output = null
    Map<String, Map> builder = null
}