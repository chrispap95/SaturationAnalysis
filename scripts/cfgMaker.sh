#!/usr/bin/sh

for i in 2000 2800
do
echo "outFilePath = out_E${i}Eta1p7.root" > simpleBH_E${i}Eta1p7.cfg
echo "filePath = root://cmseos.fnal.gov//store/user/chpapage/SingleGamma_E${i}Eta1p7/SingleGamma_E${i}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/"`ls /eos/uscms/store/user/\
chpapage/SingleGamma_E${i}Eta1p7/SingleGamma_E${i}Eta1p7_CMSSW_10_6_3_patch1_upgrade2023_D41_ntuples/`"/0000" >> simpleBH_E${i}Eta1p7.cfg
cat >> simpleBH_E${i}Eta1p7.cfg << "EOF"
recoFileName = ntuples
nRuns = 200
EOF
echo "deadfrac = 0.${df}" >> simpleBH_E${i}Eta1p7.cfg
done
