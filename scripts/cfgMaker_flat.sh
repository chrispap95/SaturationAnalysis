#!/usr/bin/sh

# Basic configuration values - most are self explanatory
#     - filesToProcess: number of files to use for each dead fraction.
#     - splitLevel: maximum number of files per job.
#                   For no splitting, set it >= ${deadFractions[0]}
energyRange=500to3000
eta=1p62
phi=0p0
pGenerator=SingleGamma
cmssw=${CMSSW_VERSION}
geometry=upgrade2023_D41
siteUrl=root://cmseos.fnal.gov/
filesToProcess=300
splitLevel=300
firstFile=0

# Prepare arrays for job splitting
configuration=${pGenerator}_E${energyRange}Eta${eta}Phi${phi}
if [[ ${phi} == Flat ]];
then
  configuration=${pGenerator}_E${energyRange}Eta${eta}
fi

numberOfJobs=0
lastJobFiles=0
if [ $((${filesToProcess}%${splitLevel})) -eq 0 ]
then
  numberOfJobs=$((${filesToProcess}/${splitLevel}))
  lastJobFiles=${splitLevel}
else
  numberOfJobs=$((${filesToProcess}/${splitLevel} + 1))
  lastJobFiles=$((${filesToProcess}%${splitLevel}))
fi

# Loop over dead fractions and jobs per dead fraction
# nRuns: number of files to process
i=0
for j in `seq ${numberOfJobs}`
do
  if [ ${numberOfJobs} -eq 1 ]
  then
    namestring=E${energyRange}Eta${eta}Phi${phi}
    if [[ ${phi} == Flat ]];
    then
      namestring=E${energyRange}Eta${eta}
    fi
  else
    namestring=E${energyRange}Eta${eta}Phi${phi}_${j}
    if [[ ${phi} == Flat ]];
    then
      namestring=E${energyRange}Eta${eta}_${j}
    fi
  fi
  samplesPath=store/user/${USER}/${configuration}/${configuration}_${cmssw}_${geometry}_ntuples/
  echo "outFilePath = out_${namestring}.root" > sampleCreator_${namestring}.cfg
  echo "filePath = ${siteUrl}/${samplesPath}"`ls /eos/uscms/${samplesPath}`"/0000" >> sampleCreator_${namestring}.cfg
  echo "recoFileName = ntuples" >> sampleCreator_${namestring}.cfg
  echo "firstRun = ${firstFile}" >> sampleCreator_${namestring}.cfg
  if [ ${j} -eq ${numberOfJobs} ]
  then
    echo "nRuns = ${lastJobFiles}" >> sampleCreator_${namestring}.cfg
    firstFile=$((${firstFile} + ${lastJobFiles}))
  else
    echo "nRuns = ${splitLevel}" >> sampleCreator_${namestring}.cfg
    firstFile=$((${firstFile} + ${splitLevel}))
  fi
  ((i+=1))
done
