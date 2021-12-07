from td_model import ThermodynamicModel
from sklearn.linear_model import LinearRegression
import numpy as np

linReg = LinearRegression()

def dict2tdm(model, treat_as = "36N"):
    a = {}
    for k in model:
        if k in ["DataIDs","Layout","en.scale"]: continue
        if isinstance(model[k], dict) and treat_as in model[k]:
            a[k] = model[k][treat_as]
        else:
            a[k] = model[k]
    a["RBSthreshold"] = a["ThDict"]
    del a["ThDict"]
    return ThermodynamicModel(a)

def evaluate_model(model_, bricks_, delta_mu, detection_th, loglums, weights=None, ax=None, c=None, forceLinear=False):
    pons_ = model_.bricks2pons({k:bricks_[k]-delta_mu for k in bricks_})
    x = np.log10(pons_)
    xmin = min(detection_th, np.percentile(x[np.isfinite(x)],.01))
    y = loglums.copy()
    if weights is None:
        weights = np.ones(len(loglums))
    if ax is not None:
        #ax.scatter(x, y, s=10*weights**.5, alpha= .5,c=c)
        ax.scatter(x, y, s=2, alpha= .5,c=c)
        #ax.axvline(detection_th, color="r", ls="--")
    x = np.maximum(x,detection_th)
    linReg.fit(x.reshape(-1,1), y, sample_weight=weights)
    if forceLinear:
        linReg.coef_[:] = 1
        linReg.intercept_ = ((y*weights).sum()-(x*weights).sum())/weights.sum()
    xr = np.linspace(detection_th,x.max())
    score = linReg.score(x.reshape(-1,1), y, sample_weight=weights)
    if ax is not None:
        ypr = linReg.predict(xr.reshape(-1,1))
        ax.plot(xr, ypr, "C3", label=r"$r=%.3f$"%score)
        ax.plot([xmin, detection_th],[ypr[0]]*2,"C3")
    return score

def find_detection_threshold(model_, bricks_, delta_mu, det_ths, loglums, weights=None, forceLinear=False):
    scores = []
    for det_th in det_ths:
        scores += [evaluate_model(model_, bricks_, delta_mu, det_th, loglums, weights=weights, forceLinear=forceLinear)]
    return {
        "detection_th_opt": det_ths[np.argmax(scores)],
        "score_opt": np.max(scores),
        "all_res": (det_ths, scores)
    }

def find_delta_mu(model_, bricks_, delta_mus, loglums, weights=None, ax=None, forceLinear=False):
    scores = []
    detths = []
    for delta_mu in delta_mus:
        pons_ = model_.bricks2pons({k:bricks_[k]-delta_mu for k in bricks_})
        logpons = np.log10(pons_)
        endpoint = np.percentile(logpons,99.5)
        det_ths = np.linspace(logpons.min(), endpoint)
        tmp = find_detection_threshold(model_, bricks_, delta_mu, det_ths, loglums, weights=weights, forceLinear=forceLinear)
        if ax is not None:
            ax.plot(*tmp["all_res"], label=delta_mu)
        scores += [evaluate_model(model_, bricks_, delta_mu, tmp["detection_th_opt"], loglums, weights=weights, forceLinear=forceLinear)]
        detths += [tmp["detection_th_opt"]]
    return {
        "delta_mu_opt": delta_mus[np.argmax(scores)],
        "detection_th_opt": detths[np.argmax(scores)],
        "score_opt": np.max(scores),
        "all_res": (delta_mus, scores)
    }