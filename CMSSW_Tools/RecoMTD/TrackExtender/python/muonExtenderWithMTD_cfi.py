import FWCore.ParameterSet.Config as cms

from RecoMTD.TrackExtender.PropagatorWithMaterialForMTD_cfi import *
from RecoMTD.TrackExtender.muonExtenderWithMTDBase_cfi import *

muonExtenderWithMTD = muonExtenderWithMTDBase.clone()
