from toil.common import Toil
from toil.job import Job

from LsJob import LsJob


if __name__ == '__main__':
    options = Job.Runner.getDefaultArgumentParser().parse_args()

    job1 = LsJob(path="/", displayName='sysFiles')
    job2 = LsJob(path="/Users/raymondleclair", displayName='userFiles')
    job3 = LsJob(path="/Users/Shared")

    job1.addChild(job2)
    job2.addChild(job3)

    with Toil(options) as toil:
        if not toil.options.restart:
            toil.start(job1)
        else:
            toil.restart()
