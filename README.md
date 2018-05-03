# SexaQuark Analysis

## Check-out code
```
mkdir ./Analysis/HexaQuark
cd ./Analysis/HexaQuark
cmsrel CMSSW_9_2_8
cd CMSSW_9_2_8/src
git clone https://github.com/gflouris/HexaAnalysis
cmsenv
scram b -j2
```

## Produce Particle Gun MC
```
cd ../..
mkdir ./ParticleGun
cd ./ParticleGun

cmsrel 
cp ../CMSSW_8_0_21/src/HexaAnalysis/TreeProducer/scripts/mc_ParticleGun/SUS-RunIISummer15GS-00146_GENSIM_cfg.py ./CMSSW_7_1_20_patch3/src/
cd CMSSW_7_1_20_patch3/src
cmsenv
cmsRun SUS-RunIISummer15GS-00146_GENSIM_cfg.py

cd ../..
cmsrel CMSSW_8_0_21
cp ../CMSSW_8_0_21/src/HexaAnalysis/TreeProducer/scripts/mc_ParticleGun/SUS-RunIISummer16DR80Premix-00068_* CMSSW_8_0_21/src/
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

## Documentation
### Adaptive Vertex Reco
- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAdaptiveVertexReconstructor
- http://cds.cern.ch/record/1166320/files/NOTE2008_033.pdf
- http://iopscience.iop.org/article/10.1088/0954-3899/34/12/N01/pdf


## MC_Analysis  --  fpartous

MonteCarlo analysis for Master thesis.

main analysis code in MCAnalysis.C
### how to run:

root -l -q -b 'MCAnalysis.C+("c1")'

where "c1" can be replaced by "c2", "c3" or "c123", depending on the category one wants to produce plots of. c1 = prompt, c2 = material interaction, c3 = decays of baryons/mesons, c123 = all three together.

outputs root file with plots

## running all three categories:

bash MCAnalysis_all_categories.bash

## exporting rootfile plots to pdf

bash exportplots_all_categories.bash


