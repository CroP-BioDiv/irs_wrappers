#!/usr/bin/env python

import sys 
from math import lgamma
from math import log

data    = open(sys.argv[1])
repeats = open(sys.argv[2])

chloro    = {'LSC' : [], 'SSC' : [] }
chlorosize =0

def lpbinom(x,n,p):
    lprob = log(p) * x + log(1-p) * (n-x) + lgamma(n+1) - lgamma(x+1) - lgamma(n-x+1)
    return lprob

# We scan the blast matches:
#    We build a vector with one position per base pair counting the matches

# The match file has the following format:
#    LSC/SSC  begin end  same_strand=1/diff_strand=0

for line in data:
    parts = line.strip().split()
    if len(parts) >= 4:
        single      = parts[0]
        begin       = int(parts[1])
        end         = int(parts[2])
        direction   = int(parts[3])
        
        # Change the code of the direction:
        #    reverse complement = -1
        if direction==0:
            direction=-1
            
        if end > chlorosize:
            extsize =  end - chlorosize 
            chloro['LSC'].extend([0] * extsize)
            chloro['SSC'].extend([0] * extsize)
            chlorosize=len(chloro['LSC'])
        
        begin-=1
        
        chr = chloro[single]
        
        for p in range(begin,end):
            chr[p]+=direction
 
maxSSC = float(max(abs(n) for n in chloro['SSC']))
maxLSC = float(max(abs(n) for n in chloro['LSC']))

chloro['SSC']=[n / maxSSC for n in chloro['SSC']]
chloro['LSC']=[n / maxLSC for n in chloro['LSC']]

scoreMax=0
len1Max=0
len2Max=0

imax = len(chloro['LSC'])

for line in repeats:
    parts   = line.strip().split()
    
    # First repeat position and length 
    # (position start at 0)
    pos1    = int(parts[1]) -1
    len1    = int(parts[3])
    
    # Second repeat position and length
    # (position start at 0)
    pos2    = int(parts[2]) -1
    len2    = int(parts[4])
    
    # Location of the central single copy 
    #       - in between the two IR -
    c_begin = min(pos1 + len1,imax)
    c_end   = min(pos2,imax)
    
    # Location of the external single copy 
    #       - in between the two IR -
    o_max   = min(pos1 ,imax)
    o_min   = min(pos2 + len2, imax)
    
    # count of coherent matches for LSC and SSC on the central single copy 
    c_lsc   = abs(sum(chloro['LSC'][n] for n in range(c_begin,c_end)))
    c_ssc   = abs(sum(chloro['SSC'][n] for n in range(c_begin,c_end)))

    # count of coherent matches for LSC and SSC on the external single copy 
    #    this score is in two parts before the first copy and after the second
    o_lsc   = sum(chloro['LSC'][n] for n in range(0,o_max))
    o_ssc   = sum(chloro['SSC'][n] for n in range(0,o_max))

    o_lsc  += sum(chloro['LSC'][n] for n in range(o_min,imax))
    o_ssc  += sum(chloro['SSC'][n] for n in range(o_min,imax))
    
    o_lsc   = abs(o_lsc)
    o_ssc   = abs(o_ssc)
    
    c = float(c_lsc + c_ssc)
    o = float(o_lsc + o_ssc)
        
    if c > 0:
        c_lsc /= c
        c_ssc /= c 
    
    if o > 0:
        o_lsc /= o 
        o_ssc /= o 
    
    score = ((c_lsc - c_ssc) ** 2 + (o_lsc - o_ssc) ** 2) / 2.0
    # pvalue= 
    # print >>sys.stderr,"c.lsc = %f c.ssc = %f   o.lsc = %f o.ssc = %f score = %6.4f (len=%d)" % (c_lsc,c_ssc,o_lsc,o_ssc,score,len1)
        
    if (score >= scoreMax) and ((len1 > len1Max) or (len2 > len2Max)):
        scoreMax = score
        pos1Max  = pos1
        pos2Max  = pos2
        len1Max  = len1
        len2Max  = len2
        
c_begin = min(pos1Max + len1Max,imax)
c_end   = min(pos2Max,imax)
o_max   = min(pos1Max,imax)
o_min   = min(pos2Max + len2Max,imax)

c_lsc   = sum(chloro['LSC'][n] for n in range(c_begin,c_end))
c_ssc   = sum(chloro['SSC'][n] for n in range(c_begin,c_end))

o_lsc   = sum(chloro['LSC'][n] for n in range(0,o_max))
o_ssc   = sum(chloro['SSC'][n] for n in range(0,o_max))

o_lsc  += sum(chloro['LSC'][n] for n in range(o_min,imax))
o_ssc  += sum(chloro['SSC'][n] for n in range(o_min,imax))

if abs(c_lsc) > abs(c_ssc):
    center = "LSC"  
    dcenter= "+" if c_lsc > 0 else "-"
else:
    center = "SSC"
    dcenter= "+" if c_ssc > 0 else "-"

if abs(o_lsc) > abs(o_ssc):
    out    = "LSC"  
    dout   = "+" if o_lsc > 0 else "-"
else:
    out    = "SSC"
    dout   = "+" if o_ssc > 0 else "-"
    
    
    
sys.stdout.write("%s %s %s %s %d %d %d %d %6.5f\n" % (center,
                                                      dcenter,
                                                      out,
                                                      dout,
                                                      pos1Max + 1,
                                                      len1Max,
                                                      pos2Max + 1,
                                                      len2Max,
                                                      scoreMax))
    
    
         
#for p in range(chlorosize):
#    sys.stdout.write("%d %d %d\n"  % (p,chloro['SSC'][p],chloro['LSC'][p]))
    
