#!/usr/bin/sh

for i in `seq 0 499`
do
echo "outFilePath = out_E0to3000Eta1p7_${i}.root" > sampleCreator_E0to3000Eta1p7_${i}.cfg
echo "filePath = root://cmseos.fnal.gov//store/user/chpapage/SingleGamma_E0to3000Eta1p7/SingleGamma_E0to3000Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/"`ls /eos/uscms/store/user/\
chpapage/SingleGamma_E0to3000Eta1p7/SingleGamma_E0to3000Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/`"/0000" >> sampleCreator_E0to3000Eta1p7_${i}.cfg
cat >> sampleCreator_E0to3000Eta1p7_${i}.cfg << "EOF"
recoFileName = ntuples
nRuns = 1
EOF
echo "firstRun = $(($i * 1))" >> sampleCreator_E0to3000Eta1p7_${i}.cfg
done
