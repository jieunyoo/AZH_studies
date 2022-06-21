# vim: set sts=4 sw=4 fdm=syntax fdl=2 et:
#
#     | ___ /   |  / _)
#     |   _ \   ' /   |  __ \
#     |    ) |  . \   |  |   |
#    _| ____/  _|\_\ _| _|  _|
#
#

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.modules.common.collectionMerger import collectionMerger
from LatinoAnalysis.NanoGardener.framework.BranchMapping import mappedOutputTree, mappedEvent

import math
from itertools import combinations, permutations

class l3KinProducer(Module):
    l3KinDefault = -9999
    Zmass = 91.1876
    Wmass = 80.4
    WH3l_btagWP = -0.5884
    newbranches = {
        'WH3l_ZVeto'     : (["F"], {}),
        'WH3l_flagOSSF'  : (["O"], {}),
        'WH3l_njet'      : (["I"], {}),
        'WH3l_nbjet'     : (["I"], {}),
        'WH3l_mtlmet'    : (["F"], {'n':3}),
        'WH3l_dphilmet'  : (["F"], {'n':3}),
        'WH3l_mOSll'     : (["F"], {'n':3}),
        'WH3l_drOSll'    : (["F"], {'n':3}),
        'WH3l_ptOSll'    : (["F"], {'n':3}),
        'WH3l_chlll'     : (["I"], {}),
        'WH3l_mlll'      : (["F"], {}),
        'WH3l_ptlll'     : (["F"], {}),
        'WH3l_ptWWW'     : (["F"], {}),
        'WH3l_mtWWW'     : (["F"], {}),
        'WH3l_dphilllmet': (["F"], {}),
        'WH3l_ptW'       : (["F"], {}),

        # for ZH3l, "l" in these variables *always* refers to the lepton not associated with the Z
        #'AZH_Wmass_Hadronic'   : (["F"], {}),
        #'AZH_Wmass_NonHadronic'   : (["F"], {}),

        'AZH_Recon_Wmass_Hadronic_Jets'  : (["F"], {}),
        'AZH_Recon_Wmass_NonHadronic' : (["F"], {}),

#note: in your histo, hadd the two hadronic ones together, and two leptonic ones
        'AZH_HadronicTopMass' : (["F"], {}),
        'AZH_HadronicTbarMass' : (["F"], {}),
        'AZH_LeptonicTopMass' : (["F"], {}),
        'AZH_LeptonicTbarMass' : (["F"], {}),


#        'AZH_GeneratorTopMass': (["F"], {}),
 #       'AZH_GeneratorTbarMass': (["F"], {}),

        #'AZH_len_CleanJetCollection': (["I"], {}),
        #'AZH_len_CleanJet_4vec' : (["I"], {}),
        #'checkZH3l_isOK' : (["I"], {}),
        #'AZH_len_bJetCollection' : (["I"], {}),
        #'AZH_len_nonbJetCollection' : (["I"], {}),
        #'AZH_isAZH' : (["I"], {}),
        #'AZH_passSelection':(["I"], {}),
        #'bJetMinValue': (["F"], {}),
        #'countTopW_Hadronic': (["I"], {}),
        #'countTBarW_Hadronic':  (["I"], {}),
        #'countTopLeptonic':  (["I"], {}),
        #'countTbarLeptonic': (["I"], {}),
    }

    def __init__(self, branch_map=''):
        self._branch_map = branch_map

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = mappedOutputTree(wrappedOutputTree, mapname=self._branch_map)

        for nameBranchKey, newBranchOpt in self.newbranches.items() :
            self.out.branch(nameBranchKey, *newBranchOpt[0], **newBranchOpt[1]);

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def _WH3l_isOk(self):
        """If a good WH3l candidate event?"""
        if not self.l3_isOk:
            return False
        sign = lambda a: 1 if a >= 0 else -1
        if abs(sum([sign(l[1]) for l in self.Lepton_4vecId])) > 1:
            return False
        return True

    def WH3l_ZVeto(self):
        """Return min mass difference in OSSF lepton pairs"""
        if not self.WH3l_isOk:
            return -1*self.l3KinDefault

        minmllDiffToZ = -1*self.l3KinDefault
        for iLep, jLep in combinations(self.Lepton_4vecId, 2):
            if iLep[1]+jLep[1] == 0:
                mllDiffToZ = abs((iLep[0]+jLep[0]).M()-self.Zmass)
                minmllDiffToZ = mllDiffToZ if mllDiffToZ < minmllDiffToZ else minmllDiffToZ
        return minmllDiffToZ

    def WH3l_flagOSSF(self):
        """Return True if OSSF lepton pair is found"""
        if not self.WH3l_isOk:
            return False

        for iLep, jLep in combinations(self.Lepton_4vecId,2):
            if iLep[1]+jLep[1] == 0:
                return True

        return False

    def WH3l_njet(self):
        if not self.WH3l_isOk:
            return self.l3KinDefault
        return sum([1 if j[0].Pt() > 40 and abs(j[0].Eta()) < 4.7 else 0 for j in self.CleanJet_4vecId])

    def WH3l_nbjet(self):
        if not self.WH3l_isOk:
            return self.l3KinDefault
        return sum([ 1 if 40 > j[0].Pt() > 20 and abs(j[0].Eta()) < 4.7 and j[1] > self.WH3l_btagWP else 0 for j in self.CleanJet_4vecId ])

    def WH3l_mtlmet(self):
        """https://en.wikipedia.org/wiki/Transverse_mass with m_lepton=0, m_met=0. """
        if not self.WH3l_isOk:
            return [self.l3KinDefault]*3
        calc_mtlmet = lambda lvec, mvec: math.sqrt(2 * lvec.Pt() * mvec.Pt() * (1 - math.cos(abs(lvec.DeltaPhi(mvec))))) if lvec.Pt() > 0 else self.l3KinDefault
        return [ calc_mtlmet(l[0], self.MET) for l in self.Lepton_4vecId ]

    def WH3l_dphilmet(self):
        """Return dphi between lepton and MET"""
        if not self.WH3l_isOk:
            return [self.l3KinDefault]*3
        return [ abs(l[0].DeltaPhi(self.MET)) for l in self.Lepton_4vecId ]

    def WH3l_mOSll(self):
        """Return mass of OS lepton pair"""
        if not self.WH3l_isOk:
            return [self.l3KinDefault]*3
        return [ (iLep[0]+jLep[0]).M() if iLep[1]*jLep[1] < 0 else self.l3KinDefault for iLep, jLep in combinations(self.Lepton_4vecId,2)]

    def WH3l_drOSll(self):
        """Return dr of OS lepton pair"""
        if not self.WH3l_isOk:
            return [self.l3KinDefault]*3
        return [ iLep[0].DeltaR(jLep[0]) if iLep[1]*jLep[1] < 0 else self.l3KinDefault for iLep, jLep in combinations(self.Lepton_4vecId,2)]

    def WH3l_ptOSll(self):
        """Return pt of OS lepton pair"""
        if not self.WH3l_isOk:
            return [self.l3KinDefault]*3
        return [ (iLep[0]+jLep[0]).Pt() if iLep[1]*jLep[1] < 0 else self.l3KinDefault for iLep, jLep in combinations(self.Lepton_4vecId,2)]

    def WH3l_chlll(self):
        """Return charge sum of leptons"""
        if not self.WH3l_isOk:
            return self.l3KinDefault
        sign = lambda a: 1 if a >= 0 else -1
        return sum([ sign(l[1]) for l in self.Lepton_4vecId])

    def WH3l_mlll(self):
        """Return invariant of leptons"""
        if not self.WH3l_isOk:
            return self.l3KinDefault
        return (self.Lepton_4vecId[0][0]+self.Lepton_4vecId[1][0]+self.Lepton_4vecId[2][0]).M()

    def WH3l_ptlll(self):
        """Return pt of leptons"""
        if not self.WH3l_isOk:
            return self.l3KinDefault
        return (self.Lepton_4vecId[0][0]+self.Lepton_4vecId[1][0]+self.Lepton_4vecId[2][0]).Pt()

    def WH3l_ptWWW(self):
        """Return pt of WH"""
        if not self.WH3l_isOk:
            return self.l3KinDefault
        return (self.Lepton_4vecId[0][0]+self.Lepton_4vecId[1][0]+self.Lepton_4vecId[2][0]+self.MET).Pt()

    def WH3l_mtWWW(self):
        """Return mt of WH"""
        if not self.WH3l_isOk:
            return self.l3KinDefault
        return math.sqrt(2*self.WH3l_ptlll()*self.MET.Pt()*(1. - math.cos(self.WH3l_dphilllmet())))

    def WH3l_dphilllmet(self):
        """Return mt of WH"""
        if not self.WH3l_isOk:
            return self.l3KinDefault
        return abs((self.Lepton_4vecId[0][0]+self.Lepton_4vecId[1][0]+self.Lepton_4vecId[2][0]).DeltaPhi(self.MET))

    def WH3l_ptW(self):
        """Return pt of lepton least likely to be from the Higgs (proxy for associated W)"""
        WH3l_ptW = self.l3KinDefault
        if self.WH3l_isOk:
            mindR = -1*self.l3KinDefault
            for iLep, jLep, kLep in permutations(self.Lepton_4vecId, 3):
                if iLep[1]*jLep[1] < 0:
                    dR = iLep[0].DeltaR(jLep[0])
                    if dR < mindR:
                        mindR = dR
                        WH3l_ptW = kLep[0].Pt()
        return WH3l_ptW

    def _ZH3l_setXLepton(self):
        """Find the lepton least likely to be part of the Z pair by invariant mass.  Double-duty as a check for OSSF pair"""
        if not self.WH3l_isOk:
            return False

        minmllDiffToZ = -1*self.l3KinDefault
        self.ZH3l_XLepton = None

        for iLep, jLep, kLep in permutations(self.Lepton_4vecId, 3):
            if iLep[1]+jLep[1] == 0:
                mllDiffToZ = abs((iLep[0]+jLep[0]).M()-self.Zmass)
                if mllDiffToZ < minmllDiffToZ:
                    self.ZH3l_XLepton = kLep
                    self.pTZ = (iLep[0]+jLep[0]).Pt()
                    self.checkZmass = (iLep[0]+jLep[0]).M()
                    minmllDiffToZ = mllDiffToZ
        # print self.ZH3l_XLepton
        # print (self.ZH3l_XLepton is not None)
        return self.ZH3l_XLepton is not None


    ##code added below - jieun****************************************************##
    def returnDeltaRValue(self, quarkCandidate):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            minimumDeltaR = 9999
            self.closestCleanJetCandidate = ROOT.TLorentzVector()
            for j in self.AZH_nonbJet_GetWMass:
                 delR_computation = quarkCandidate.DeltaR(j)
                 if delR_computation < minimumDeltaR:
                    minimumDeltaR = delR_computation
                    self.closestCleanJetCandidate = j
        return minimumDeltaR

    def returnDeltaRIndex(self, quarkCandidate):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            minimumDeltaR = 9999
            matchedIndex = 0
            self.closestCleanJetCandidate = ROOT.TLorentzVector()
            #for j in self.AZH_nonbJet_GetWMass:
            for j,k in enumerate(self.AZH_nonbJet_GetWMass):
                 delR_computation = quarkCandidate.DeltaR(k)
                 if delR_computation < minimumDeltaR:
                    minimumDeltaR = delR_computation
                    self.closestCleanJetCandidate = k
                    matchedIndex = j
        return matchedIndex

    def computeDeltaR(self, quarkCandidate):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            minimumDeltaR = 9999
            self.closestCleanJetCandidate = ROOT.TLorentzVector()
            for j in self.AZH_nonbJet_GetWMass:
                 delR_computation = quarkCandidate.DeltaR(j)
                 if delR_computation < minimumDeltaR:
                    minimumDeltaR = delR_computation
                    self.closestCleanJetCandidate = j
            if minimumDeltaR > 0.4:
                 self.deltaRcutCleanJet = False
                 return self.l3KinDefault
            return self.closestCleanJetCandidate
        else:
            pass

    def AZH_Recon_Wmass_Hadronic_Jets(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.deltaRcutCleanJet is False:
            return self.l3KinDefault
        if self.IsAZH is True:
            jvec0 = self.closestCleanJetCandidate1
            jvec1 = self.closestCleanJetCandidate2
            return (jvec0 + jvec1).M()
        else:
            return self.l3KinDefault

    def AZH_Recon_Wmass_NonHadronic(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            jvec0 = self.ZH3l_XLepton[0]
            jvec1 = self.reconNeutrino
            return (jvec0 + jvec1).M()
        else:
            return self.l3KinDefault

    def computeDeltaR_Bjet(self, quarkCandidate):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            minimumDeltaR = 9999
            self.closestCleanJetCandidate = ROOT.TLorentzVector()
            for j in self.AZH_bJetList:
                 delR_computation = quarkCandidate.DeltaR(j)
                 if delR_computation < minimumDeltaR:
                    minimumDeltaR = delR_computation
                    self.closestCleanJetCandidate = j
            return minimumDeltaR

    #separate definition for b jets because 2 diff. sets of  that jets that can fail deltaR cut
    def computeDeltaR_Truth_TopB_Jet(self, quarkCandidate):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            minimumDeltaR = 9999
            self.closestBJetCandidate = ROOT.TLorentzVector()
            for j in self.AZH_bJetList:
                 delR_computation = quarkCandidate.DeltaR(j)
                 if delR_computation < minimumDeltaR:
                    minimumDeltaR = delR_computation
                    self.closestBJetCandidate = j
            if minimumDeltaR > 0.4:
                 self.deltaRcutTopBQuark = False
                 return self.l3KinDefault
            return self.closestBJetCandidate
        else:
            pass

    def computeDeltaR_Truth_TbarB_Jet(self, quarkCandidate):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            minimumDeltaR = 9999
            self.closestBJetCandidate = ROOT.TLorentzVector()
            for j in self.AZH_bJetList:
                 delR_computation = quarkCandidate.DeltaR(j)
                 if delR_computation < minimumDeltaR:
                    minimumDeltaR = delR_computation
                    self.closestBJetCandidate = j
            if minimumDeltaR > 0.4:
                 self.deltaRcutTbarBQuark = False
                 return self.l3KinDefault
            return self.closestBJetCandidate
        else:
            pass

#graphs for MC fits ********************************************************************************
#need to HADD AZH_HadronicTopMass and AZH_TbarMass to get "Hadronic Top" and similarly for leptonic top
    def AZH_HadronicTopMass(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.deltaRcutTopBQuark is False:
            return self.l3KinDefault
        if self.deltaRcutCleanJet is False:
            return self.l3KinDefault
        if self.IsAZH is True:
            if self.quarksFromWPlus is True:
                jvec0 = self.closestCleanJetCandidate1
                jvec1 = self.closestCleanJetCandidate2
                jvec2 = self.bquark1ClosestJet
                return (jvec0 + jvec1 + jvec2).M()
            else:
                return self.l3KinDefault
        else:
            return self.l3KinDefault

    def AZH_HadronicTbarMass(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.deltaRcutTbarBQuark is False:
            return self.l3KinDefault
        if self.IsAZH is True:
            if self.countquarksFromWminus is True:
                jvec0 = self.closestCleanJetCandidate1
                jvec1 = self.closestCleanJetCandidate2
                jvec2 = self.bquark2ClosestJet
                return (jvec0 + jvec1 + jvec2).M()
            else:
                jvec0 = self.ZH3l_XLepton[0]
                jvec1 = self.MET
                jvec2 = self.bquark2ClosestJet
                return (jvec0 + jvec1 + jvec2).M()
        else:
            return self.l3KinDefault

    def AZH_LeptonicTopMass(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.deltaRcutTopBQuark is False:
            return self.l3KinDefault
        if self.deltaRcutCleanJet is False:
            return self.l3KinDefault
        if self.IsAZH is True:
            if self.TopLeptons is True:
                jvec0 = self.ZH3l_XLepton[0]
                jvec1 = self.reconNeutrino
                jvec2 = self.bquark1ClosestJet
                return (jvec0 + jvec1 + jvec2).M()
            else:
                return self.l3KinDefault
        else:
            return self.l3KinDefault

    def AZH_LeptonicTbarMass(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.deltaRcutTbarBQuark is False:
            return self.l3KinDefault
        if self.deltaRcutCleanJet is False:
            return self.l3KinDefault
        if self.IsAZH is True:
            if self.TbarLeptons is True:
                jvec0 = self.ZH3l_XLepton[0]
                jvec1 = self.reconNeutrino
                jvec2 = self.bquark2ClosestJet
                return (jvec0 + jvec1 + jvec2).M()
            else:
                return self.l3KinDefault
        else:
            return self.l3KinDefault

#***********************************************************************************************************

#************below uses actual generator level quarks for comparison purpsoes
    def AZH_Wmass_Hadronic(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            jvec0 = self.quark1
            jvec1 = self.quark2
            return (jvec0 + jvec1).M()
        else:
            return self.l3KinDefault

    def AZH_Wmass_NonHadronic(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            return (self.neutrino + self.lepton).M()
        else:
            return self.l3KinDefault

    def AZH_GeneratorTopMass(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            return self.topMass
        else:
            return self.l3KinDefault

    def AZH_GeneratorTbarMass(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH is True:
            return self.tbarMass
        else:
            return self.l3KinDefault

#cutflow checks
    def checkZH3l_isOK(self):
        if self.ZH3l_isOk:
            return 1
        else:
            return 0

    def AZH_len_CleanJet_4vec(self):
        if not self.ZH3l_isOk:
            return self.l3KinDefault
        return len(self.AZH_CleanJet_4vecId)

    def AZH_len_bJetCollection(self):
        if not self.ZH3l_isOk:
            return self.l3KinDefault
        if not self.pass4JetCut:
            return self.l3KinDefault
        return len(self.bJetCollection)

    def AZH_len_nonbJetCollection(self):
        if not self.ZH3l_isOk:
            return self.l3KinDefault
        if not self.pass4JetCut:
            return self.l3KinDefault
        return len(self.AZH_nonbJet)

    def AZH_len_CleanJetCollection(self):
        if not self.ZH3l_isOk:
            return self.l3KinDefault
        return len(self.ZH3l_CleanJet_4vecId)

    def AZH_passSelection(self):
        if self.passSelection:
            return 1
        else:
            return 0

    def AZH_isAZH(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH:
            return 1
        else:
            return 0

    def bJetMinValue(self):
        if not self.passSelection:
            return self.l3KinDefault
        if self.IsAZH:
            return self.deltaRBjetMin
        else:
            return self.l3KinDefault

    def countTopW_Hadronic(self):
        if not self.passSelection:
            return self.l3KinDefault
        if not self.IsAZH:
            return self.l3KinDefault
        if self.quarksFromWPlus is True:
            return 1
        else:
            return self.l3KinDefault

    def countTBarW_Hadronic(self):
        if not self.passSelection:
            return self.l3KinDefault
        if not self.IsAZH:
            return self.l3KinDefault
        if self.quarksFromWminus is True:
            return 1
        else:
            return self.l3KinDefault

    def countTopLeptonic(self):
        if not self.passSelection:
            return self.l3KinDefault
        if not self.IsAZH:
            return self.l3KinDefault
        if self.TopLeptons is True:
            return 1
        else:
            return self.l3KinDefault

    def countTbarLeptonic(self):
        if not self.passSelection:
            return self.l3KinDefault
        if not self.IsAZH:
            return self.l3KinDefault
        if self.TbarLeptons is True:
            return 1
        else:
            return self.l3KinDefault


#******************************************************************************************
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        event = mappedEvent(event, mapname=self._branch_map)
        # Order in pt the collection merging muons and electrons
        # lepMerger must be already called
        Lepton = Collection(event, "Lepton")
        self.Lepton_4vecId = []
        for iLep in range(3):
            if len(Lepton) > iLep:
                self.Lepton_4vecId.append((ROOT.TLorentzVector(), Lepton[iLep].pdgId))
                self.Lepton_4vecId[-1][0].SetPtEtaPhiM(Lepton[iLep].pt, Lepton[iLep].eta, Lepton[iLep].phi, 0)

        self.MET = ROOT.TLorentzVector()
        self.MET.SetPtEtaPhiM(event.PuppiMET_pt, 0, event.PuppiMET_phi, 0)

        self.CleanJet_4vecId = []
        Jet = Collection(event, "Jet")
        for j in Collection(event, "CleanJet"):
            #self.CleanJet_4vecId.append((ROOT.TLorentzVector(), Jet[j.jetIdx].btagCMVA))
            self.CleanJet_4vecId.append((ROOT.TLorentzVector(), Jet[j.jetIdx].btagDeepB))
            self.CleanJet_4vecId[-1][0].SetPtEtaPhiM(j.pt, j.eta, j.phi, 0)

            print(j.jetIdx)


        self.l3_isOk = False if len(self.Lepton_4vecId) < 3 else True
        self.WH3l_isOk = self._WH3l_isOk()
        self.ZH3l_isOk = self._ZH3l_setXLepton()

##code added below - jieun****************************************************
        self.ZH3l_CleanJet_4vecId = []
        self.AZH_CleanJet_4vecId = []
        self.bJetCollection = []
        self.AZH_nonbJet = []
        self.pass4JetCut = False
        self.passBJetCut = False
        self.passSelection = False

        if self.ZH3l_isOk:
            self.ZH3l_CleanJet_4vecId = [ j for j in self.CleanJet_4vecId if j[0].Pt() > 30 and abs(j[0].Eta()) < 4.7]
            self.AZH_CleanJet_4vecId = [ j for j in self.ZH3l_CleanJet_4vecId if len(self.ZH3l_CleanJet_4vecId) >= 4]
            if len(self.AZH_CleanJet_4vecId) >= 4: self.pass4JetCut = True
            if self.pass4JetCut == True:
                self.bJetCollection = [ j for j in self.AZH_CleanJet_4vecId if j[1] > 0.4941]
                self.AZH_nonbJet = [j for j in self.AZH_CleanJet_4vecId if j[1] <= 0.4941]

                if len(self.bJetCollection) >=2:
                    self.passBJetCut = True

                #these separate lists are to make it easier to grab LorentzVector
                self.AZH_nonbJet_GetWMass = []
                for i in self.AZH_nonbJet:
                    self.AZH_nonbJet_GetWMass.append(i[0])
                self.AZH_bJetList = []
                for i in self.bJetCollection:
                    self.AZH_bJetList.append(i[0])

        if self.ZH3l_isOk and self.pass4JetCut and self.passBJetCut:
            self.passSelection = True
        else:
            self.passSelection = False

        #get true quarks for matching
        genPartList = Collection(event, "GenPart")
        self.genPart_list = []
        self.countquarks = 0
        self.countleptons = 0
        self.IsAZH = False
        self.deltaRcutCleanJet = True
        self.deltaRcutTopBQuark = True
        self.deltaRcutTbarBQuark = True
        self.countquarksFromWplus = 0
        self.quarksFromWPlus = False
        self.countquarksFromWminus = 0
        self.quarksFromWminus = False
        self.countTopLeptons = 0
        self.countTbarLeptons = 0
        self.TopLeptons = False
        self.TbarLeptons = False

        if self.passSelection:
            for j in genPartList:
                if ((j.statusFlags>>7 &1) and j.genPartIdxMother >= 0 and abs(genPartList[j.genPartIdxMother].pdgId) == 24):
                    if abs(j.pdgId) <= 4:
                       self.countquarks = self.countquarks+1
                    if (abs(j.pdgId) == 13 or abs(j.pdgId) == 11 or abs(j.pdgId) == 12 or abs(j.pdgId) == 14 or abs(j.pdgId) == 15 or abs(j.pdgId) == 16 and abs(genPartList[j.genPartIdxMother].pdgId) == 24) :
                       self.countleptons = self.countleptons+1
                if ((j.statusFlags>>7 &1) and j.genPartIdxMother >= 0 and abs(genPartList[j.genPartIdxMother].pdgId) == 6):
                    if abs(j.pdgId) == 5:
                        self.countquarks = self.countquarks+1

            if (self.countquarks == 4 and self.countleptons == 2):
                self.IsAZH = True
            else:
                self.IsAZH = False

        twoQuarkList_eta = []
        twoQuarkList_phi = []
        twoQuarkList_pt = []
        twoQuarkList_mass = []
        neutrinoEta = 0
        neutrinoPt = 0
        neutrinoPhi = 0
        leptonEta = 0
        leptonPt = 0
        leptonMass = 0
        leptonPhi = 0
        topMass = 0
        topPt = 0
        topPhi = 0
        topEta = 0
        tbarMass = 0
        tbarPt =  0
        tbarPhi = 0
        tbarEta = 0
        self.topMass = 0
        self.tbarMass = 0
        bquark1Mass = 0 #this is the b-quark
        bquark1Pt = 0
        bquark1Phi = 0
        bquark1Eta = 0
        bquark2Mass = 0 #this is bbar
        bquark2Pt = 0
        bquark2Phi = 0
        bquark2Eta = 0

        if self.IsAZH is True and self.passSelection:
            for j in genPartList: #to get the W mass
                if ( (j.statusFlags>>7 &1) and j.genPartIdxMother >= 0 and abs(genPartList[j.genPartIdxMother].pdgId) == 24):
                    if ( (abs(j.pdgId) == 12 or abs(j.pdgId) == 14 or abs(j.pdgId) == 16 and j.status == 1 )):
                        neutrinoPt = j.pt
                        neutrinoEta = j.eta
                        neutrinoPhi = j.phi
                    if ( (abs(j.pdgId) == 11 or abs(j.pdgId) == 13 or abs(j.pdgId) == 15 and j.status == 1)):
                        leptonPt = j.pt
                        leptonEta = j.eta
                        leptonPhi = j.phi
                        leptonMass = j.mass
                    if ( abs(j.pdgId) < 5):
                        twoQuarkList_eta.append(j.eta)
                        twoQuarkList_phi.append(j.phi)
                        twoQuarkList_pt.append(j.pt)
                        twoQuarkList_mass.append(j.mass)

                    if(genPartList[j.genPartIdxMother].pdgId) == 24: #check to see if mother is a W+, if so the hadronic quarks go with W+
                        if ( abs(j.pdgId) < 5):
                            self.countquarksFromWplus += 1
                        if ( (abs(j.pdgId) == 11 or abs(j.pdgId) == 13 or abs(j.pdgId) == 15 and j.status == 1)):
                            self.countTopLeptons +=1
                        if ( (abs(j.pdgId) == 12 or abs(j.pdgId) == 14 or abs(j.pdgId) == 16 and j.status == 1 )):
                            self.countTopLeptons +=1
                    if(genPartList[j.genPartIdxMother].pdgId) == -24: #check to see if mother is a W-
                        if ( abs(j.pdgId) < 5):
                            self.countquarksFromWminus += 1
                        if ( (abs(j.pdgId) == 11 or abs(j.pdgId) == 13 or abs(j.pdgId) == 15 and j.status == 1)):
                            self.countTbarLeptons +=1
                        if ( (abs(j.pdgId) == 12 or abs(j.pdgId) == 14 or abs(j.pdgId) == 16 and j.status == 1 )):
                            self.countTbarLeptons += 1

                #this is just to check the generator value of the top mass
                if ( (j.statusFlags>>7 &1) and j.genPartIdxMother >= 0 and abs(genPartList[j.genPartIdxMother].pdgId) == 35):
                   if ( j.pdgId == 6):
                       topMass = j.mass
                       #print(topMass)
                       topPt = j.pt
                       topPhi = j.phi
                       topEta = j.eta
                   if ( j.pdgId == -6):
                       tbarMass = j.mass
                       tbarPt = j.pt
                       tbarPhi = j.phi
                       tbarEta = j.eta
                if ( (j.statusFlags>>7 &1) and j.genPartIdxMother >= 0 and genPartList[j.genPartIdxMother].pdgId == 6):
                    #bottom quark has pdgID == 5
                    if ( j.pdgId == 5):
                        bquark1Mass = j.mass
                        bquark1Pt = j.pt
                        bquark1Phi = j.phi
                        bquark1Eta = j.eta
                if ( (j.statusFlags>>7 &1) and j.genPartIdxMother >= 0 and genPartList[j.genPartIdxMother].pdgId == -6):
                    if ( j.pdgId == -5):
                        bquark2Mass = j.mass
                        bquark2Pt = j.pt
                        bquark2Phi = j.phi
                        bquark1Eta = j.eta

            if self.countquarksFromWplus == 2:
                self.quarksFromWPlus = True
            if self.countquarksFromWminus  == 2:
                self.quarksFromWminus = True
            if self.countTopLeptons == 2:
                self.TopLeptons = True
            if self.countTbarLeptons == 2:
                self.TbarLeptons = True

            #quark1 and quark2 are for hadronic jets
            self.quark1 = ROOT.TLorentzVector()
            self.quark2 = ROOT.TLorentzVector()
            self.neutrino = ROOT.TLorentzVector()
            self.lepton = ROOT.TLorentzVector()
            self.top = ROOT.TLorentzVector()
            self.tBar = ROOT.TLorentzVector()
            self.topMass = topMass
            self.tbarMass = tbarMass
            self.bquark1 = ROOT.TLorentzVector()
            self.bquark2 = ROOT.TLorentzVector()
            self.deltaRquark1 = 0
            self.deltaRquark2 = 0
            self.deltaRbquark1 = 0
            self.deltaRbquark2 = 0

            self.quark1.SetPtEtaPhiM(twoQuarkList_pt[0], twoQuarkList_eta[0], twoQuarkList_phi[0], twoQuarkList_mass[0])
            self.quark2.SetPtEtaPhiM(twoQuarkList_pt[1], twoQuarkList_eta[1], twoQuarkList_phi[1], twoQuarkList_mass[1])
            self.neutrino.SetPtEtaPhiM(neutrinoPt, neutrinoEta, neutrinoPhi, 0.)
            self.lepton.SetPtEtaPhiM(leptonPt, leptonEta, leptonPhi, leptonMass)
            self.top.SetPtEtaPhiM(topPt,topEta,topPhi,topMass)
            self.tBar.SetPtEtaPhiM(tbarPt,tbarEta,tbarPhi,tbarMass)
            self.bquark1.SetPtEtaPhiM(bquark1Pt,bquark1Eta,bquark1Phi,bquark1Mass)
            self.bquark2.SetPtEtaPhiM(bquark2Pt,bquark2Eta,bquark2Phi,bquark2Mass)

            #match the quarks to the closest cleanJet
            self.closestCleanJetCandidate1 = ROOT.TLorentzVector()
            self.closestCleanJetCandidate2 = ROOT.TLorentzVector()
            self.closestCleanJetCandidate1 = self.computeDeltaR(self.quark1)
            self.closestCleanJetCandidate2 = self.computeDeltaR(self.quark2)




            #match the bquarks to the closest clean jets - bquark1 is from top, bquark2 is from tbar
            self.bquark1ClosestJet = ROOT.TLorentzVector() #bquark1 = b
            self.bquark2ClosestJet = ROOT.TLorentzVector() #bquark2 = bbar
            self.bquark1ClosestJet = self.computeDeltaR_Truth_TopB_Jet(self.bquark1)
            self.bquark2ClosestJet = self.computeDeltaR_Truth_TbarB_Jet(self.bquark2)

            #check fractions
            self.checkDR1 = 0
            self.checkDR2 = 0
            self.indexDR1 = 0
            self.indexDR2 = 0
            self.deltaRBjetMin = 0
            self.neutrinopz1 = 0
            self.neutrinopz2 = 0
            self.reconNeutrino = ROOT.TLorentzVector()

            if self.IsAZH is True and self.passSelection:
                  self.checkDR1 = self.returnDeltaRValue(self.quark1)
                  self.indexDR1 = self.returnDeltaRIndex(self.quark1)
                  self.checkDR2 = self.returnDeltaRValue(self.quark2)
                  self.indexDR2 = self.returnDeltaRIndex(self.quark2)
                  self.deltaRBjetMin = self.computeDeltaR_Bjet(self.bquark1)

                  #TO DO: change this into a function but for now this is fine
                  #neutrino recon
                  if (self.ZH3l_XLepton):
                    zeta = 0.5 * pow(80.4,2) + self.MET.Pt()*self.ZH3l_XLepton[0].Pt()*math.cos(self.ZH3l_XLepton[0].DeltaPhi(self.MET))
                    mess = (pow(zeta,2)*pow(self.ZH3l_XLepton[0].Pz(),2)) / pow(self.ZH3l_XLepton[0].Pt(),4)  -  (pow(self.MET.Pt(),2)*pow(self.ZH3l_XLepton[0].E(),2) - pow(zeta,2)) / pow(self.ZH3l_XLepton[0].Pt(),2)
                    if mess > 0:
                        mess = math.sqrt(mess)
                    else:
                        mess = 0
                    self.neutrinopz1 = (zeta * self.ZH3l_XLepton[0].Pz())/pow(self.ZH3l_XLepton[0].Pt(),2) + mess
                    self.neutrinopz2 = (zeta * self.ZH3l_XLepton[0].Pz())/pow(self.ZH3l_XLepton[0].Pt(),2) - mess

                    self.minNeutrino = 9999
                    self.chooseThisNeutrinoPZ = 9999

                    self.minDiff1 = self.neutrinopz1 - self.neutrino.Pz()
                    self.minDiff2 = self.neutrinopz2 - self.neutrino.Pz()
                    if abs(self.minDiff1) < abs(self.minDiff2):
                        self.minNeutrino = self.minDiff1
                        self.chooseThisNeutrinoPZ = self.neutrinopz1
                    else:
                        self.minNeutrino = self.minDiff2
                        self.chooseThisNeutrinoPZ = self.neutrinopz2

                    self.neutrinoEnergy = math.sqrt( pow(self.MET.Pt(),2) + pow(self.chooseThisNeutrinoPZ,2) )
                    #print('true, guess1, guess2', self.neutrino.Pz(), self.neutrinopz1, self.neutrinopz2)
                    self.reconNeutrino.SetPxPyPzE(self.MET.Px(),self.MET.Py(),self.chooseThisNeutrinoPZ,self.neutrinoEnergy)


        for nameBranchKey in self.newbranches.keys():
            self.out.fillBranch(nameBranchKey, getattr(self, nameBranchKey)());

        return True
