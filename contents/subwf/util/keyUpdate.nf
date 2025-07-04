include {keyJoin} from './keyJoin.nf'

workflow keyUpdate {
        take:
        left
        right
        options

        main:
        options = options instanceof Map ? options + [update: true] : [by: options, update: true]
        keyJoin(left, right, options) | set{result}

        emit:
        result
}