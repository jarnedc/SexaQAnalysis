# HexaQuark Analysis

## Check-out code
```
mkdir ./Analysis/HexaQuark
cd ./Analysis/HexaQuark
cmsrel CMSSW_8_0_21
cd CMSSW_8_0_21/src
git cms-init
git clone https://github.com/gflouris/HexaAnalysis
scram b -j2
```

## Produce Particle Gun MC
```
cd ../..
mkdir ./ParticleGun
cd ./ParticleGun

cmsrel CMSSW_7_1_20_patch3
cd CMSSW_7_1_20_patch3/src
cmsenv
cmsRun SUS-RunIISummer15GS-00146_GENSIM_cfg.py

cd ../..
cmsrel CMSSW_8_0_21
cmsenv
cd CMSSW_8_0_21/src
cmsenv
cmsRun SUS-RunIISummer16DR80Premix-00068_Step1_cfg.py
cmsRun SUS-RunIISummer16DR80Premix-00068_Step2_cfg.py
cd ../../..
```

## Produce ntuples
```
cd ./Analysis/HexaQuark/CMSSW_8_0_21/src/TreeProducer/test
cmsenv
cmsRun treeproducer_AOD_MC_cfg.py
```

## Template macro
```
cd ../macros/
root -l analysis.C
```
