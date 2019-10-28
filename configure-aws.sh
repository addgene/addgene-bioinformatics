sudo yum update -y
sudo yum install gcc -y

# install Docker
sudo amazon-linux-extras install docker
sudo service docker start
sudo usermod -a -G docker ec2-user

# install Conda's Python 3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
sudo ./Miniconda3-latest-Linux-x86_64.sh -b -p /usr/miniconda
export PATH="/usr/miniconda/bin/:$PATH"
conda update -y -n base conda && conda config --add channels conda-forge && conda config --add channels defaults &&  conda config --add channels bioconda

# install Toil itself
pip install toil --user
