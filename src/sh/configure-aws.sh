# Install essentials
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
        emacs           \
        wget            \
        zip             \
        build-essential \
        pkg-config      \
        bzip2           \
        ca-certificates

# Install Docker (from https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-16-04)
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
sudo apt-get update
apt-cache policy docker-ce
sudo apt-get install -y docker-ce
sudo usermod -a -G docker $USER

# Install APT's pip and virtualenv
sudo apt-get install python-pip
sudo apt-get install python-virtualenv
sudo apt-get install virtualenvwrapper
