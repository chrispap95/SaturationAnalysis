#!/usr/bin/sh

# Basic configuration values - self explanatory

pGenerator=CloseBySingleGamma
cmssw=${CMSSW_VERSION}
geometry=upgrade2026_D54
siteUrl=root://cmseos.fnal.gov/
eta=1p62
phi=0p0
energies=(10 50 100 200 300 400 500 600 700 900 1100 1300 1600 2000 2400 2900)
filesToProcess=100

# Loop over dead fractions and energies
# nRuns: number of files to process

for En in ${energies[@]}
do
  if [ ${phi} -eq Flat ];
  then
    namestring=E${En}Eta${eta}
  else
    namestring=E${En}Eta${eta}Phi${phi}
  fi
  configuration=${pGenerator}_${namestring}
  samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/
  echo "outFilePath = out_${namestring}.root" > sampleCreator_${namestring}.cfg
  echo "filePath = ${siteUrl}/${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> sampleCreator_${namestring}.cfg
  echo "recoFileName = ntuples" >> sampleCreator_${namestring}.cfg
  echo "nRuns = ${filesToProcess}" >> sampleCreator_${namestring}.cfg
done
