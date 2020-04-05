#!/usr/bin/sh

for i in `seq 0 19`
do
for E in 400 550 750 1000 1400
do
echo "outFilePath = out_E${E}Eta1p7_${i}.root" > sampleCreator_E${E}Eta1p7_${i}.cfg
echo "filePath = root://cmseos.fnal.gov//store/user/chpapage/SingleGamma_E${E}Eta1p7/SingleGamma_E${E}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/"`ls /eos/uscms/store/user/\
chpapage/SingleGamma_E${E}Eta1p7/SingleGamma_E${E}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/`"/0000" >> sampleCreator_E${E}Eta1p7_${i}.cfg
cat >> sampleCreator_E${E}Eta1p7_${i}.cfg << "EOF"
recoFileName = ntuples
nRuns = 4
EOF
echo "firstRun = $(($i * 4))" >> sampleCreator_E${E}Eta1p7_${i}.cfg
done
done
