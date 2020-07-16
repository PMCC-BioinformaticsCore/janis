# Requirements for using Janis

Janis has a number of requirements for building and running workflows.

## tL;dr - Simple

- Unix machine (macOS / Ubuntu / RHEL / CentOS / etc)
- Docker
- Python 3.6+
- Zip (system unarchive utility)

This is the easiest configuration without any extra software.

## Simple requirements

To just build pipelines (without running them), you don't need Docker or Singularity.

Most engines (including CWLTool / Cromwell) work with

## Software description 

### CWLTool

CWLTool has the following requirements:


## Setting up an AWS EC2 to run Janis:


```bash
sudo yum update -y
sudo yum install -y yum-utils
sudo amazon-linux-extras install -y docker
sudo service docker start
sudo usermod -a -G docker ec2-user
sudo yum install -y python3
sudo yum install -y gcc python3-devel

python3 -m venv ~/janis/env
source ~/janis/env/bin/activate
pip install janis-pipelines

# reconnect to get new permissions
exit
```

