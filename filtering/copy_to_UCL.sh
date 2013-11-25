# Script for copying a directory from homedir at ULB to storage element at UCL. You need to be able to connect to ingrid via e.g. ssh ingrid-ui1.cism.ucl.ac.be. It is recommended to establish password-less connection as detailed in https://cp3.irmp.ucl.ac.be/projects/cp3admin/wiki/UsersPage/Support/CP3/IT/SshKeys or in https://jez.web.cern.ch/jez/SshKeys.html (if you do not access to CP3 Twiki)

# Example usage: source copy_to_UCL.sh ZHttNtuples/53X/Data/DoubleElectron_Run2012D-PromptReco-v1/
# The rsync commands checks the content of the target directory at UCL and synchronises it with the local directory, so even if copy command fails for various reasons (timeout...), it can be painlessly resumed

rsync -av ${1}/ ingrid-ui1.cism.ucl.ac.be:/storage_rw/data/cms/store/user/${USER}/${1}/