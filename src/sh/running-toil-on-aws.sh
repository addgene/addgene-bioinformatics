
# TODO: Try this again
# See:
#   https://toil.readthedocs.io/en/latest/running/cloud/amazon.html#prepareaws
#   https://toil.readthedocs.io/en/latest/gettingStarted/quickStart.html#launching-a-toil-workflow-in-aws

# aws configure
# emacs ~/.boto

TOIL_APPLIANCE_SELF=quay.io/ucsc_cgl/toil:3.18.0

toil launch-cluster clustername \
     --zone us-east-1a \
     --keyPairName id_rsa \
     --leaderNodeType t2.medium

toil rsync-cluster clustername \
     --zone us-east-1a \
     ../python/helloWorld.py :/tmp

toil ssh-cluster clustername \
     --zone us-east-1a 

# python /tmp/helloWorld.py aws:us-east-1:helloWorld-s3-bucket

toil destroy-cluster clustername \
     --zone us-east-1a 
