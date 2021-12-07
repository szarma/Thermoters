import numpy as np
from general_functions import multi_map
try: from scipy.special import logsumexp
except: from scipy.misc import logsumexp
from fastFunctions import tensum, bindingEnergies, getDiNu

from model_functions import getBricks
class ThermodynamicModel:
    def __init__(self, parameters):
        self.params = dict(parameters)

    def sequences2bricks(self, seqs, dinuCoordsAndValues=None):
        dinucl = dinuCoordsAndValues is not None
        if dinucl:
            dinuCoords, dinuValues = dinuCoordsAndValues
        strands = ["frw"]
        if self.params["includeRC"]:
            strands += ["rc"]
        bricks = {}
        for strand in strands:
            sq = seqs
            if strand=="rc":
                sq = 3 - sq[:, ::-1].copy(order="C")
            tmp = getBricks(
                self.params["matrices"],
                self.params["min.spacer"],
                self.params["sp.penalties"],
                sq,
                makeLengthConsistent=True).T
            tmp += -self.params["chem.pot"]

            if dinucl:
                global mp_getDiNu

                def mp_getDiNu(coord_):
                    return getDiNu(*coord_,
                                   n1=self.params["matrices"][0].shape[0],
                                   minSpacer=self.params["min.spacer"],
                                   n2=self.params["matrices"][1].shape[0],
                                   sequences=sq,
                                   nSpacer=len(self.params["sp.penalties"])).T

                tmpDn = np.array(multi_map(mp_getDiNu, dinuCoords, processes=14))
                #                 tmpDn = np.array([
                #                     getDiNu(*coord,
                #                         n1=mdl["matrices"][0].shape[0],
                #                         minSpacer=mdl["min.spacer"],
                #                         n2=mdl["matrices"][1].shape[0],
                #                         sequences=sq,
                #                         nSpacer=len(mdl["sp.penalties"])).T
                #                     for coord in dinuCoords])
                tmp += np.array(tensum(dinuValues, tmpDn))
            if strand=="rc":
                tmp = tmp[:, ::-1]
            bricks[strand] = tmp
        return bricks

    def bricks2pons(self,bricks):
        if "logClearanceRate" in self.params:
            R_ = np.exp(self.params["logClearanceRate"])
        else:
            R_ = 0
        bdni = bricks["frw"]
        thresholdPos = self.params["RBSthreshold"]
        if thresholdPos <= 0:
            thresholdPos = bdni.shape[1] + thresholdPos
        off = thresholdPos < bdni.shape[1]
        bindMode_ = self.params["bindMode"]
        if bindMode_ == "add":
            bindsumF = lambda xi: np.sum(
                1. / (np.exp(xi) + R_),
                axis=tuple(range(1, xi.ndim))
            )
        if bindMode_ == "max":
            bindsumF = lambda xi: np.exp(-np.min(xi, axis=tuple(range(1, xi.ndim))))
        sumON_ = bindsumF(bdni[:, :thresholdPos])
        if off:
            sumOFF_ = bindsumF(bdni[:, thresholdPos:])
        else:
            sumOFF_ = 0.
        if "rc" in bricks:
            bdni = bricks["rc"]
            rcOcclusion = self.params.get("rcOcclusion", np.arange(bdni.shape[1]))
            sumOFF_ += bindsumF(bdni[:, rcOcclusion])
        Pons_ = sumON_ / (1. + sumOFF_ + sumON_)
        return Pons_
    def __repr__(self):
        return self.params.__repr__()