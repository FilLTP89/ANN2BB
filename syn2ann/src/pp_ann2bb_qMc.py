# -*- coding: utf-8 -*-
#!/usr/bin/env python
u'''Plot tools for nn'''
u'''Required modules'''
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import pandas as pd
import seaborn as sns
import argparse
from os import listdir
from os.path import isfile, join
import numpy as np


u'''General informations'''
__author__ = "Filippo Gatti & Didier Clouteau"
__copyright__ = "Copyright 2018, CentraleSupélec (MSSMat UMR CNRS 8579)"
__credits__ = ["Filippo Gatti"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Filippo Gatti"
__email__ = "filippo.gatti@centralesupelec.fr"
__status__ = "Beta"
def setup():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p','--pwd',default='/tmp1/gattif/ann_sobol',help='path to dataset')
    parser.add_argument('-c','--cpt',default='all',nargs = '*',help='station labels')
    parser.add_argument('-d','--drs',default='all',nargs = '*',help='direction')
    opt = parser.parse_args()
    return opt
    
if __name__=='__main__':

    opt = setup()
    opt.cpt = map(str.lower,opt.cpt)
    opt.drs = map(str.lower,opt.drs)
    print(opt.cpt)
    if 'all' in opt.cpt:
        if 'all' in opt.drs:
            fnm = [f for f in listdir(opt.pwd) if isfile(join(opt.pwd, f))]
        else:
            fnm = [f for f in listdir(opt.pwd) if isfile(join(opt.pwd, f)) and any(d in f.lower() for d in opt.drs)]
    else:
        if 'all' in opt.drs:
            fnm = [f for f in listdir(opt.pwd) if isfile(join(opt.pwd, f)) and any(c in f.lower() for c in opt.cpt)]
        else:
            fnm = [f for f in listdir(opt.pwd) if isfile(join(opt.pwd, f)) and any(c in f.lower() for c in opt.cpt) and any(d in f.lower() for d in opt.drs)]
    nfn = len(fnm)
    idx = [int(f.split('.txt')[0].split('_')[1]) for f in fnm]
    fnm = [fnm[idx.index(i)] for i in range(0)
    print(fnm)

    exit()
    xf=[]
    df=[]
    af=[]
    for i in fnm:
        if 'X_' in i:
            xf.append(i)
        elif 'AA_' in i:
            af.append(i)
        elif 'D_' in i:
            df.append(i)
    nr = []
    for i in range(10,10+len(xf)):
        xn=xf[xf.index('X_{:d}.txt'.format(i), )]
        dn=df[df.index('D_{:d}.txt'.format(i), )]
        an=af[af.index('AA_{:d}.txt'.format(i), )]
        xv=np.loadtxt(join(pwd,xn), np.float_,',')
        dv=np.loadtxt(join(pwd,dn), np.float_,',')
        av=np.loadtxt(join(pwd,an), np.float_,',')
        av=av.reshape((int(np.sqrt(av.size)),int(np.sqrt(av.size))))
        xe = np.linalg.solve(av,dv)
        nr.append(np.sum(xe**2-xv**2))
        print(nr[-1])
    exit()

