#!/usr/bin/sh
source ${PWD}/prepareCondor.sh

eta=1p7
energies=(5 10 20 40 60 80 100)

for En in ${energies[@]}
do
namestring=E${En}Eta${eta}
argument=sampleCreator_${namestring}.cfg\ out_${namestring}.root\ ${CMSSW_VERSION}\ ${USER}

# Write jdl file
cat > condor_${namestring}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
EOF
echo "Transfer_Input_Files = condor-exec.csh, ${CMSSW_VERSION}.tgz" >> condor_${namestring}.jdl
echo "Arguments = ${argument}" >> condor_${namestring}.jdl
cat >> condor_${namestring}.jdl << "EOF"
Output = sampleCreator_$(Cluster)_$(Process).stdout
Error = sampleCreator_$(Cluster)_$(Process).stderr
Log = sampleCreator_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF

# Submit job
condor_submit condor_${namestring}.jdl
done
