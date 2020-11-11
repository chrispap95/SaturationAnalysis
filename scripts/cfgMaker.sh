#!/usr/bin/sh

# Basic configuration values - self explanatory

eta=1p7
pGenerator=SingleGamma
cmssw=${CMSSW_VERSION}
geometry=upgrade2023_D41
siteUrl=root://cmseos.fnal.gov/
energies=(5 10 20 40 60 80 100)
filesToProcess=200

# Loop over dead fractions and energies
# nRuns: number of files to process

for En in ${energies[@]}
do
  namestring=E${En}Eta${eta}
  configuration=${pGenerator}_E${En}Eta${eta}
  samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/
  echo "outFilePath = out_${namestring}.root" > sampleCreator_${namestring}.cfg
  echo "filePath = ${siteUrl}/${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> sampleCreator_${namestring}.cfg
  echo "recoFileName = ntuples" >> sampleCreator_${namestring}.cfg
  echo "nRuns = ${filesToProcess}" >> sampleCreator_${namestring}.cfg
done
