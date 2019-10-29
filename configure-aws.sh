sudo apt-get clean all && \
    sudo apt-get update && \
    sudo apt-get upgrade -y && \
    sudo apt-get install -y  \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        build-essential \
        pkg-config      \
        bzip2           \
        ca-certificates && \
    sudo apt-get clean && \
    sudo apt-get purge && \
    sudo rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install Docker (from https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04)
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt-get update
apt-cache policy docker-ce
sudo apt-get install -y docker-ceq

# install Conda's Python 3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sudo bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/bin
miniconda3/bin/conda init bash
source ~/.bashrc

# install Toil itself
pip install toil --user
