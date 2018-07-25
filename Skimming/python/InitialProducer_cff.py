import FWCore.ParameterSet.Config as cms

InitialProducer = cms.EDProducer(
  'InitialProducer',
  TrackCollection       = cms.InputTag("generalTracks"),
)

