## Definitions:
from collections import OrderedDict
import numpy as np
# from copy import deepcopy
# from general_functions import tally,mode
from general_functions import multi_map

try:
    from scipy.special import logsumexp
except:
    from scipy.misc import logsumexp

from fastFunctions import tensum, bindingEnergies, getDiNu
    
def slideSingleMatrix(m, seqs):
    Lout = seqs.shape[1]-m.shape[0]+1
    return np.array([bindingEnergies(m,seqs[:,offset:offset+m.shape[0]]) for offset in range(Lout)]).T

def getBricks(twoMatrices,
              minSpacer,
              spacerPenalties,
              sequences,
              makeLengthConsistent=False):

    n1, n2 = [m.shape[0] for m in twoMatrices]
    nSpacer = len(spacerPenalties)
    nSeq, seqL = sequences.shape
    
    energyBoxes = np.array([
        slideSingleMatrix(twoMatrices[0], sequences[:,              : seqL-n2-minSpacer]).T,
        slideSingleMatrix(twoMatrices[1], sequences[:, n1+minSpacer :                  ]).T
                  ])
    if makeLengthConsistent:    
        # works only if the center spacer Penalty is 0
        spFlex = nSpacer//2
        assert spacerPenalties[spFlex]==0
        Lmatrix = minSpacer+n1+n2+spFlex
        Lbrick = seqL-Lmatrix+1
        effergies = np.ones((nSpacer,Lbrick,nSeq))*100  # large number so it vanishes when exp(-#)
        for iS in range(nSpacer):
            try:
                tmp =  energyBoxes[0,:energyBoxes.shape[1]-iS]\
                                   +energyBoxes[1,iS:]\
                                   +spacerPenalties[iS]
    #             print (tmp.shape,Lbrick)
                tmp = tmp[-Lbrick:]
                effergies[iS][-tmp.shape[0]:] = tmp
            except:
                pass
    else:
        effergies = [np.array(  energyBoxes[0,:energyBoxes.shape[1]-iS] \
                               +energyBoxes[1,iS:]
                               +spacerPenalties[iS])
                     for iS in range(nSpacer)]
        # to align the sequences by right-flushing
        effergies = np.array([effergies[iS][nSpacer-iS:] for iS in range(nSpacer)])

    return effergies

def getBrickDict(seqDict,mdl,dinucl=False,
                 subtractChemPot=True,
                 useChemPot="chem.pot",
                 makeLengthConsistent=False,
                 dinuCoordsAndValues = None):
    if dinucl:
        dinuCoords, dinuValues = dinuCoordsAndValues
    out = OrderedDict()
    strands = [0]
    if mdl["includeRC"]:
        strands += [1]
    for did in seqDict:
        for strand in strands:
            sq = seqDict[did]
            if strand:
                sq = 3-sq[:,::-1].copy(order="C")
            tmp = getBricks(
                mdl["matrices"],
                mdl["min.spacer"],
                mdl["sp.penalties"],
                sq,
                makeLengthConsistent=makeLengthConsistent).T
            if subtractChemPot:
                try:
                    mu = mdl[useChemPot][did]
                except:
                    for k in mdl[useChemPot]:
                        if did in k:
                            mu = mdl[useChemPot][k]
                            break
                tmp += -mu
                
            if dinucl:
                global mp_getDiNu
                def mp_getDiNu(coord_):
                    return getDiNu(*coord_,
                        n1=mdl["matrices"][0].shape[0],
                        minSpacer=mdl["min.spacer"],
                        n2=mdl["matrices"][1].shape[0],
                        sequences=sq,
                        nSpacer=len(mdl["sp.penalties"])).T
                    
                tmpDn = np.array(multi_map(mp_getDiNu, dinuCoords, processes = 14))
#                 tmpDn = np.array([
#                     getDiNu(*coord,
#                         n1=mdl["matrices"][0].shape[0],
#                         minSpacer=mdl["min.spacer"],
#                         n2=mdl["matrices"][1].shape[0],
#                         sequences=sq,
#                         nSpacer=len(mdl["sp.penalties"])).T
#                     for coord in dinuCoords])
                tmp += np.array(tensum(dinuValues,tmpDn))
            if strand:
                tmp = tmp[:,::-1]
            out[did+"_rc"*strand] = tmp
    return out


    
def brick2lps(bricks_DNIs,
              fitpars,
              thresholdPosDict_ = None,
              bindMode_ = None,
              useChemPot = "chem.pot"
             ):
    out = {}               
    if thresholdPosDict_ is None:
        thresholdPosDict_ = fitpars["ThDict"]
    if bindMode_ is None:
        bindMode_ = fitpars["bindMode"]
    try:
        R_ = np.exp(fitpars["logClearanceRate"])
    except:
        R_ = None
    for dataID_ in bricks_DNIs:
        if "_rc" in dataID_: continue
        
        bdni = bricks_DNIs[dataID_]
        
        try:
            thresholdPos = thresholdPosDict_.get(dataID_,thresholdPosDict_["Prl"])
        except:
            for k in thresholdPosDict_:
                if dataID_ in k:
                    thresholdPos = thresholdPosDict_[k]
                    break
        if thresholdPos<=0:
            thresholdPos = bdni.shape[1]+thresholdPos

        off = thresholdPos<bdni.shape[1]
        if bindMode_ == "add":
            if R_ is None:
                bindF = lambda xi: -logsumexp(-xi,axis=tuple(range(1,xi.ndim)))
            else:
                bindF = lambda xi: -np.log(np.sum(
                                    1./(np.exp(xi)+R_),
                                    axis=tuple(range(1,xi.ndim))
                                            ))
        if bindMode_ == "max":
            bindF = lambda xi: np.min(xi,axis=tuple(range(1,xi.ndim)))
        effON_  =     bindF(bdni[:,:thresholdPos])
        if off:
            effOFF_ = bindF(bdni[:,thresholdPos:])
        else:
            effOFF_ = 0.
        if dataID_+"_rc" in bricks_DNIs:
            bdni_rc = bricks_DNIs[dataID_+"_rc"]
            rcOcclusion = fitpars.get("rcOcclusion",np.arange(bdni_rc.shape[1]))
            effOFF_ += bindF(bdni_rc[:,rcOcclusion])
        Pons_ = np.exp(-effON_)/(1.+np.exp(-effON_)+np.exp(-effOFF_))
        out[dataID_] = np.log10(Pons_)
    return out
    


def lps2eval(fitpar, objF, numData,
             DataIDs_   = None,
             tt         = "training",
             fit        = False,
             logPonDict_= None,
             bricks_    = None,
             binEdges_  = None,
             dinucl = False,
             dinuCoordsAndValues = None,
             useChemPot = None
             ):

    if DataIDs_ is None:
        DataIDs_ = fitpar["DataIDs"]
    if fit is None:
        fit = "train" in tt
    if useChemPot is None:
        try:
            fitpar["chem.pot_%s"%objF]
            useChemPot = "chem.pot_%s"%objF
        except:
            useChemPot = "chem.pot"
#         print ("using", useChemPot)
            
    data_ = numData[tt]
    if logPonDict_ is None:
        if bricks_ is None:
            bricks_ = getBrickDict( {did: data_[did]["seqs"] for did in DataIDs_}, fitpar, dinucl=dinucl, dinuCoordsAndValues = dinuCoordsAndValues, useChemPot=useChemPot)
        esc = fitpar["en.scale"]
        logPonDict_ = brick2lps( {el: bricks_[el]*esc for el in bricks_ }, fitpar)
    out = {}
        
    for dataID_ in DataIDs_:
        seqs_     =       data_[dataID_]["seqs"]
        digiLums_ =       data_[dataID_]["digiLums"]
        weights_  =       data_[dataID_]["weights"]
        Ndata_    =   len(data_[dataID_]["weights"])
        logPon_   = logPonDict_[dataID_].reshape(-1,1)
#         if not np.all(np.isfinite(logPon_)):
#             if objF=="mlogL":
#                 out[dataID_] = np.inf
#             else:
#                 out[dataID_] = np.nan
#             continue

        logisticRegression = fitpar["logisticRegression"]
        try:
            LR_ = logisticRegression[dataID_]
        except:
            for k in logisticRegression:
                if dataID_ in k:
                    LR_ = logisticRegression[k]
                    break
                    
        if fit: LR_.fit( logPon_, digiLums_, sample_weight=weights_)
            
        if objF=="mlogL":
            out[dataID_] = -(LR_.predict_log_proba(logPon_)[range(Ndata_),digiLums_]*weights_).sum()
            
        if objF=="linR2":
            from sklearn.linear_model import LinearRegression
            LinReg = LinearRegression()
            lums = data_[dataID_]["lums"]
            LinReg.fit(logPon_,lums,sample_weight=weights_)
            out[dataID_] = LinReg.score(logPon_,lums,sample_weight=weights_)

        if objF=="r2":
            lums = data_[dataID_]["lums"]
            errs = LR_.predict(logPon_)-lums
            wmse = sum(weights_*errs**2)/weights_.sum()
            lums_mean = np.sum(lums*weights_)/weights_.sum()
            wvar = sum(weights_*(lums-lums_mean)**2)/weights_.sum()
            out[dataID_] = 1-wmse/wvar

    return out 

    
def reprBigM(theFitPars):
    Lspacer = theFitPars["Layout"][1]
    ms = [m-np.repeat(m.min(axis=1),4).reshape(-1,4) for m in theFitPars["matrices"]]
    bigM = np.vstack([ms[0]]+
              [np.ones((1,4))*np.nan]*Lspacer+
              [ms[1]]
             )
    return bigM

# def showdf(a_):
#     from IPython.display import display
#     display(a_.applymap("{0:0.1f}".format).style.set_properties(**{'text-align': 'right'}))