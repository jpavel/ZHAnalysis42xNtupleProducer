voms-proxy-init --voms cms:/cms/becms --valid 100:0
cp $X509_USER_PROXY  myProxy
voms-proxy-info --all
