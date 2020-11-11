#!/usr/bin/sh
source ${PWD}/prepareCondor.sh

# Define submission parameters
#   - samplesNumber is the number of files to process per dead fraction
#     to find the proper number you need to weight the samples with the
#     dead fractions.
energyRange=0to3000
eta=1p7
numberOfJobs=2

for j in `seq ${numberOfJobs}`
do
if [ ${numberOfJobs} -eq 1 ]
then
  namestring=E${energyRange}Eta${eta}
else
  namestring=E${energyRange}Eta${eta}_${j}
fi
argument=sampleCreator_${namestring}.cfg\ out_${namestring}.root\ ${CMSSW_VERSION}\ ${USER}

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
condor_submit condor_${namestring}.jdl
done
