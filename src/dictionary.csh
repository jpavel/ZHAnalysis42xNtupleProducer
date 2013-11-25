#!/bin/csh

#eval `scramv1 runtime -csh`
rootcint -f eventdict.cc -c -I${PWD}/../../.. \
         #-p Analysis/NtupleProducer/interface/mytrack.h \
            Analysis/NtupleProducer/interface/myobject.h \
            Analysis/NtupleProducer/interface/myevent.h \
            Analysis/NtupleProducer/interface/LinkDef.h
