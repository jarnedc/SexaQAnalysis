This directory contains a few scripts:
-plotTSalis.py: to plot pt according to TSalis function, normalised
-plotPzTSalis.py: to plot a pz according to TSalis and uniform eta, normalised
-simulation.cc: c++ macro to simulate some kinematics: a particle (such as the S) hitting another particle (such as a neutron). The combined particle decays to 2 other particles (such as Ks and Lambda) and in the rest frame of this particle the 2 particles are send out isotropically and back to back. Then the LorentzVectors of these decay products are boosted to the LAB frame. Like this we have some info on which delta_phi and delta_theta distributions to expect.
