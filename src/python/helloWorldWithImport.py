from toil.common import Toil
from toil.job import Job


class HelloWorld(Job):
    def __init__(self, inputFileID):
        Job.__init__(self,  memory="2G", cores=2, disk="3G")
        self.inputFileID = inputFileID

    def run(self, fileStore):
        with fileStore.readGlobalFileStream(self.inputFileID) as fi:
            with fileStore.writeGlobalFileStream() as (fo, outputFileID):
                fo.write("==1== " + fileStore.localTempDir + '\n')
                fo.write("==2== " + fileStore.getLocalTempDir() + '\n')
                fo.write("==3== " + fileStore.localTempDir + '\n')
                fo.write("==4== " + fileStore.getLocalTempDir() + '\n')
                fo.write("==5== " + fileStore.localTempDir + '\n')
                fo.write("==6== " + fi.read() + 'World!')
            return outputFileID


if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()

    with Toil(options) as toil:
        if not toil.options.restart:

            read_one_file_id = toil.importFile(
                "file:///Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/dat/miscellaneous/A11967A_sW0154_FASTQ/A11967A_sW0154_A01_R1_001.fastq.gz",
                sharedFileName="A11967A_sW0154_A01_R1_001.fastq.gz")

            inputFileID = toil.importFile('file:///Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/src/python/prefix.txt')
            outputFileID = toil.start(HelloWorld(inputFileID))
        else:
            outputFileID = toil.restart()

        toil.exportFile(outputFileID, 'file:///Users/raymondleclair/Projects/Addgene/addgene-bioinformatics/src/python/greeting.txt')
