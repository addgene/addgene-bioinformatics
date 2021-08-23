
path = "/Users/williamspear/projects/springbok/addgene-bioinformatics/dat/simple-data/*.fastq.gz"
read_channel = Channel.fromFilePairs('/Users/williamspear/projects/springbok/addgene-bioinformatics/dat/simple-data/*_R{1,2}_*.fastq.gz')

adapters_channel = Channel.fromPath('/Users/williamspear/projects/springbok/addgene-bioinformatics/src/python/jobs/adapters.fa')

// format fq for bb processing
fq_ch = read_channel.map { arr -> [(arr[0] - ~/\w+(?<=_)/), (arr[0] - ~/_\w+/), arr[1]].flatten() }

process bbmerge {

echo true
stageInMode = 'copy'

input:
tuple val(pair_id), val(plate_id), file(read1), file(read2) from fq_ch
file(adapters) from adapters_channel

output:
path("*")

beforeScript 'docker login -u ralatsdio -p "$RALATSDIO_TOKEN"'
afterScript 'docker logout'

script:
"""
bbduk.sh > ${pair_id}.txt
"""
}

//out_ch.view {print "$it"}

//files_for_spades_assembler_ch.subscribe {println it}
