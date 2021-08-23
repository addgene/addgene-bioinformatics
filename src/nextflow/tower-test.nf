process bbmerge {

beforeScript 'docker login -u ralatsdio -p "$RALATSDIO_TOKEN"'
afterScript 'docker logout'

script:
"""
bbduk.sh > bbduk.txt
"""
}

//out_ch.view {print "$it"}

//files_for_spades_assembler_ch.subscribe {println it}
