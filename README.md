# SaturationAnalyzer
Code that performs simpleBH analysis using CMSSW produced ntuples as input.

## Introduction
Works on ```CMSSW_11_X_Y```. Setup an environment and clone into ```$CMSSW_BASE/src```.
Then, do
```bash
mkdir {bin,lib,obj}
source setup.sh
make
```

## Execute code
To execute the code
```bash
./bin/sampleCreator -c scripts/sampleCreator.cfg
```

You can view the output
```bash
root -l out.root
```

## Use condor to submit
*__Beware__: These instructions and scripts currently work only on the LPC cluster and contain my username. In order to use them, you have to put your username and if you are not on the LPC cluster you need to edit the code accordingly.*

I strongly suggest that you setup a fresh CMSSW release in order to make Condor submissions. So, do
```bash
mkdir condorSubmissions
cd condorSubmissions
cmsrel CMSSW_11_2_0_pre8
cd CMSSW_11_2_0_pre8/src
cmsenv
git clone https://github.com/chrispap95/SaturationAnalysis.git
scram b
```
Then, prepare your CMSSW to trasfer it to the nodes by issuing
```bash
cd SaturationAnalyzer
```

Then, submit a job
```bash
mkdir logs
sh condorSubmitter.sh
```
