ansible-playbook workstation.yml


Notes:

1. Make sure every remote host is configured with your SSH key so that passwordless SSH can happen for Ansible.

The inventory file contains the username that is used.

Example: ssh-copy-id bhuyan.m@vc7hpc-085.hpc.wehi.edu.au