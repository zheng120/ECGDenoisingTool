# ***************************************************************************
# Copyright 2017-2019, Jianwei Zheng, Chapman University,
# zheng120@mail.chapman.edu
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# 	http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Written by Jianwei Zheng.
from scipy import signal
import matplotlib.pyplot as plt
def NLM_1dDarbon(signal,Nvar,P,PatchHW):
    if isinstance(P,int): # scalar has been entered; expand into patch sample index vector
        P = P-1 #Python start index from 0
        Pvec = np.array(range(-P,P+1))
    else:
        Pvec = P # use the vector that has been input
    signal = np.array(signal)
    #debug = [];
    N = len(signal)

    denoisedSig = np.empty(len(signal)) #NaN * ones(size(signal));
    denoisedSig[:] = np.nan
    # to simpify, don't bother denoising edges
    iStart = PatchHW+1
    iEnd = N - PatchHW
    denoisedSig[iStart: iEnd] = 0

    #debug.iStart = iStart;
    #debug.iEnd = iEnd;

    # initialize weight normalization
    Z = np.zeros(len(signal))
    cnt = np.zeros(len(signal))

    # convert lambda value to  'h', denominator, as in original Buades papers
    Npatch = 2 * PatchHW + 1
    h = 2 * Npatch * Nvar**2

    for idx in Pvec: # loop over all possible differences: s - t
        # do summation over p - Eq.3 in Darbon
        k = np.array(range(N))
        kplus = k + idx
        igood = np.where((kplus >=0) & (kplus < N)) # ignore OOB data; we could also handle it
        SSD = np.zeros(len(k))
        SSD[igood] = (signal[k[igood]] - signal[kplus[igood]])**2
        Sdx = np.cumsum(SSD)

        for ii in range(iStart,iEnd): # loop over all points 's'
            distance = Sdx[ii + PatchHW] - Sdx[ii - PatchHW-1] #Eq 4;this is in place of point - by - point MSE
            # but note the - 1; we want to icnlude the point ii - iPatchHW

            w = math.exp(-distance/h) # Eq 2 in Darbon
            t = ii + idx # in the papers, this is not made explicit

            if t>0 and t<N:
                denoisedSig[ii] = denoisedSig[ii] + w * signal[t]
                Z[ii] = Z[ii] + w
                #cnt[ii] = cnt[ii] + 1
                #print('ii',ii)
                #print('t',t)
                #print('w',w)
                #print('denoisedSig[ii]', denoisedSig[ii])
                #print('Z[ii]',Z[ii])
     # loop over shifts

    # now apply normalization
    denoisedSig = denoisedSig/(Z + sys.float_info.epsilon)
    denoisedSig[0: PatchHW+1] =signal[0: PatchHW+1]
    denoisedSig[ - PatchHW: ] =signal[- PatchHW: ]
    #debug.Z = Z;

    return denoisedSig#,debug
