# author zouhu
import astropy.stats as ast
import numpy as np


def rebin_image(a, shape, fun=np.sum):
    ashape = a.shape
    if a.shape[0] != shape[0] or a.shape[1] != shape[1]:
        sh = shape[0], a.shape[0] // shape[0], shape[1], a.shape[1] / shape[1]
        return fun(fun(a.reshape(sh), -1), 1)
    else:
        return a


def sigmaclip_limitsig(data, error=None, sig_limit=None, **kwd):
    data = np.array(data)
    mdata = ast.sigma_clip(data, **kwd)
    if sig_limit is not None:
        while (True):
            med = np.ma.median(mdata)
            sig = np.ma.std(mdata)
            if sig < sig_limit: break
            index = np.ma.argmax(np.ma.abs(mdata - med))
            mdata.mask[index] = True
    return mdata


def bin_stats(x, y=None, minbin=None, maxbin=None, nbin=30, bins=None, **kwd):
    if minbin is None: minbin = np.min(x)
    if maxbin is None: maxbin = np.max(x)
    if nbin is None: nbin = 30
    if bins is None: bins = np.linspace(minbin, maxbin, nbin + 1)
    if y is None: y = x
    nbin = len(bins)
    yy = np.zeros((5, nbin - 1))
    for i in range(nbin - 1):
        yy[0, i] = (bins[i] + bins[i + 1]) / 2.0
        mask = (x >= bins[i]) & (x < bins[i + 1])
        yy[4, i] = mask.sum()
        if mask.any():
            ymask = sigmaclip_limitsig(y[mask], **kwd)
            yy[1, i] = np.ma.median(ymask)
            yy[2, i] = np.ma.std(ymask)
            yy[3, i] = ymask.mask.sum()
    return yy


def binnum_stats(x, y, minbin=None, maxbin=None, binnum=None, binedges=None, **kwd):
    x = np.array(x)
    y = np.array(y)
    nx = len(x)
    if len(y) != nx: raise ValueError("y should have same size of x")
    index = np.argsort(x)
    x = x[index]
    y = y[index]
    if minbin is None: minbin = np.min(x)
    if maxbin is None: maxbin = np.max(x)
    mask = (x >= minbin) & (x <= maxbin)
    x = x[mask];
    y = y[mask]
    nx = len(x)
    if binnum is None: binnum = nx // 10
    binedges = [x[0]]
    yy = []
    count = 0
    while (True):
        count += 1
        ind0 = (count - 1) * binnum
        if ind0 >= nx: break
        ind1 = count * binnum
        if ind1 > nx: ind1 = nx
        binedges.append(x[ind1 - 1])
        ibin = (binedges[count] + binedges[count - 1]) / 2.0
        slc = slice(ind0, ind1)
        ymask = sigmaclip_limitsig(y[slc], **kwd)
        yy.append([ibin, np.ma.median(ymask), np.ma.std(ymask), (~ymask.mask).sum()])
    return np.array(yy)


def int2inf(data):
    # get integer towards positive/negative infinite
    data = np.array(data)
    intdata = np.floor(data)
    mask = data > 0
    intdata[mask] = np.ceil(data[mask])
    return intdata


def valid_coordinates(coordinates, size=None):
    """
    convert tuple, list or np.ndarray of cooridinates to a 2d ndarray
      and check the coordinates within an image size
    INPUT
        coordinates: 2d array for coordinates
        size: size of an image
    OUTPUT:
        coord: reshape coordinates
        indcoord: index of coordinates in a size range
    """
    if isinstance(coordinates, (list, tuple, np.ndarray)):
        coord = np.atleast_2d(coordinates)
        if coord.shape[0] != 2 and coord.shape[1] != 2:
            raise ValueError("coordinates should have at least one axis with 2 elements")
        if coord.shape[1] != 2 and coord.shape[0] == 2:
            coord = coord.transpose()
    else:
        raise TypeError("coordinates should be list or array of (x,y) pixel positions")

    if size is None:
        return coord
    else:
        if len(size) != 2:
            raise ValueError("size should have 2 elements")
        nx, ny = size
        x = coord[:, 0]
        y = coord[:, 1]
        indcoord = np.arange(coord.shape[0])
        good = (x >= 0.5) & (x < nx + 0.5) & (y >= 0.5) & (y < ny + 0.5)
        if np.any(good):
            indcoord = indcoord[good]
        else:
            raise ValueError('coordinates are not in the image range')
        return coord, indcoord


def closest_match(coord1, coord2, min_dist=1.0):
    """
    find closest pairs between two sets of coordinates
    coord1: coordinates to be matched
    coord2: coordinates matched to
    min_dist: separation tolerance

    output: idx1,idx2
        idx1: matched index for coord1
        idx2: matched index for coord2
    """
    coord1 = valid_coordinates(coord1)
    coord2 = valid_coordinates(coord2)

    n1 = len(coord1)
    n2 = len(coord2)
    index1 = []
    index2 = []
    x2 = coord2[:, 0]
    y2 = coord2[:, 1]
    for i in range(n1):
        ix1, iy1 = coord1[i]
        index = np.where((np.abs(x2 - ix1) < min_dist) & (np.abs(y2 - iy1) < min_dist))[0]
        nmatch = len(index)
        if nmatch < 1:
            continue
        if nmatch > 1:
            x2tmp = x2[index];
            y2tmp = y2[index]
            dist = ((x2tmp - ix1) ** 2 + (y2tmp - iy1) ** 2)
            indsort = np.argsort(dist)
            index1.append(i)
            index2.append(index[indsort])
        else:
            index1.append(i)
            index2.append(index)
    return index1, index2


def weighted_mean(x, sigmax, weight_square=True):
    """
    Calculate the mean and estimated errors for a set of data points
    DESCRIPTION:
       This routine is adapted from Program 5-1, XFIT, from "Data Reduction
       and Error Analysis for the Physical Sciences", p. 76, by Philip R.
       Bevington, McGraw Hill.  This routine computes the weighted mean using
       Instrumental weights (w=1/sigma^2).
    INPUTS:
      x      - Array of data points
      sigmax - array of standard deviations for data points
      weight_square - if True, weight is invariance, else the reciprocal of the error
    OUTPUTS:
      xmean  - weighted mean
      sigmam - standard deviation of mean
      stdm - standard deviation of data
    """
    x = np.atleast_1d(x).copy()
    sigmax = np.atleast_1d(sigmax).copy()
    if len(x) == 1:
        xmean = x
        sigmam = sigmax
        stdm = sigmax
    else:
        weight = 1.0 / sigmax ** 2
        weight1 = weight
        if not weight_square: weight1 = 1.0 / sigmax
        wsum = weight.sum()
        xmean = (weight1 * x).sum() / weight1.sum()
        sigmam = np.sqrt(1.0 / wsum)
        stdm = np.sqrt(np.sum((x - xmean) ** 2) / len(x))
    return xmean, sigmam, stdm
