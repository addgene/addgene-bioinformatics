
path = "/Users/williamspear/projects/springbok/addgene-bioinformatics/dat/simple-data/*.fastq.gz"
read_channel = Channel.fromFilePairs('/Users/williamspear/projects/springbok/addgene-bioinformatics/dat/test-data/*_R{1,2}_*.fastq.gz')

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

publishDir "$workflow.projectDir/ngs-${curr_date}/${run_name}/final_fq/${plate_id}_CLEAN/", mode: "copy", pattern: "*.clean.fastq.qz"
publishDir "$workflow.projectDir/ngs-${curr_date}/${run_name}/final_fq/${plate_id}_NORM/", mode: "copy", pattern: "*.norm.fq"
publishDir "$workflow.projectDir/ngs-${curr_date}/${run_name}/final_fq/${plate_id}_MERGED/", mode: "copy", pattern: "*.merged.fq"

script:
"""
/bbmap/bbduk.sh > ${pair_id}.txt

/bbmap/bbduk.sh \
in=$read1 \
in2=$read2 \
out=${pair_id}_1.clean.fastq.gz \
out2=${pair_id}_2.clean.fastq.gz \
ref=/bbmap/resources/nextera.fa.gz \
ktrimright=t k=27 hdist=1 edist=0 qtrim=rl trimq=25 minlength=30 trimbyoverlap=t minoverlap=24 ordered=t qin=33 \
overwrite=true

/bbmap/bbnorm.sh -Xmx1G \
in=${pair_id}_1.clean.fastq.gz \
in2=${pair_id}_2.clean.fastq.gz \
out=${pair_id}_1.norm.fq \
out2=${pair_id}_2.norm.fq \
target=100 maxdepth=6 qin=33 \
overwrite=true \

/bbmap/bbmerge.sh -Xmx1G \
in=${pair_id}_1.norm.fq \
in2=${pair_id}_2.norm.fq \
out=${pair_id}.merged.fq \
outu=${pair_id}_1.unmerged.fq \
outu2=${pair_id}_2.unmerged.fq \
vloose=t qin=33 \
"""
}
