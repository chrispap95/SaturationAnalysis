#!/usr/bin/sh

for i in `seq 0 499`
do
cat > condor_E0to3000Eta1p7_${i}.jdl << "EOF"
universe = vanilla
Executable = condor-exec.csh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = condor-exec.csh, CMSSW_10_6_3_patch1.tgz
EOF
echo "Arguments = sampleCreator_E0to3000Eta1p7_${i}.cfg out_E0to3000Eta1p7_${i}.root" >> condor_E0to3000Eta1p7_${i}.jdl
cat >> condor_E0to3000Eta1p7_${i}.jdl << "EOF"
Output = sampleCreator_$(Cluster)_$(Process).stdout
Error = sampleCreator_$(Cluster)_$(Process).stderr
Log = sampleCreator_$(Cluster)_$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue 1
EOF
condor_submit condor_E0to3000Eta1p7_${i}.jdl
done
