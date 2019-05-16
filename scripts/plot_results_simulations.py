#sample to generate fisher p-values

import math
import numpy as np
import scipy.stats as stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-n1f", type=int ,help="if n1 is fixed or not",default=1)
parser.add_argument("-p1", type=float ,help="p1 prob",default=0.5)
parser.add_argument("-n1", type=float ,help="n1 prob",default=0.5)
parser.add_argument("-n", type=float ,help="n ",default=10000.)
args = parser.parse_args()

p0 = args.p1
p1 = args.p1
pn1 = args.n1
n = args.n
n1 = int(n*pn1)
n0 = int(n - n1)
pvalues_fisher = list()
pvalues_barnard = list()
fixn1 = args.n1f
#pvaluesb = list()
#pvaluesb = list()
fisher_file_path = "fisher_sampled_pvalue_"+str(n)+"_"+str(pn1)+"_"+str(p1)+".csv"
barnard_file_path = "barnard_sampled_pvalue_"+str(n)+"_"+str(pn1)+"_"+str(p1)+".csv"


fin = open(fisher_file_path,'r')
count = 0
for line in fin:
    try:
        #print line
        items = line.split(';')
        if len(items) > 2:
            if float(items[1]) == n and int(items[2]) == fixn1 and float(items[0]) > 0.0:
                pvalues_fisher.append( float(items[0]) )
                count = count + 1
    except ValueError:
        pass
print count
pvalues_fisher.sort()


fin = open(barnard_file_path,'r')
count = 0
for line in fin:
    try:
        items = line.split(';')
        if len(items) > 2:
            if float(items[1]) == n and int(items[2]) == fixn1 and float(items[0]) > 0.0:
                pvalues_barnard.append( float(items[0]) )
                count = count + 1
    except ValueError:
        pass
print count
pvalues_barnard.sort()


import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['figure.dpi'] = 900
import matplotlib.pyplot as plt
f, ax1 = plt.subplots(figsize=(6,6))
ax1.set_yscale('log', nonposy="clip")
ax1.set_xscale('log', nonposy="clip")

pointsf = np.linspace(0, 1, len(pvalues_fisher), endpoint=False)
pointsb = np.linspace(0, 1, len(pvalues_barnard), endpoint=False)
minvalue = min(pvalues_fisher[0],pvalues_barnard[0])*0.9
print minvalue
pointsline = np.linspace(minvalue, 1, len(pvalues_fisher), endpoint=False)
#print points
plt.plot(pointsline , pointsline , '-' , color='black')
plt.ylabel(r'$p$-value')
plt.xlabel(r'expected $p$-value')
ordmagn = math.floor(math.log10(float(len(pvalues_fisher))) + 0.5)
plt.title(r'$n='+str(int(args.n))+'$ , $p_{n_1}='+str(args.n1)+'$ , $p_1='+str(args.p1)+'$ , $j_p=10^'+str(int(ordmagn))+'$')
plt.scatter(pointsf , pvalues_fisher,marker='x' , s=3 , c='red', label='Fisher')
plt.scatter(pointsb , pvalues_barnard,marker='x' , s=3 , c='blue', label='Barnard')
ax1.set_xlim(1.0 , minvalue)
ax1.set_ylim(1.0 , minvalue)
plt.legend(loc='upper left')
#print pvaluesb
#plt.plot(points , pvaluesb)
plt.savefig('fisher_barnard_sim'+str(fixn1)+'_'+str(n)+'_'+str(pn1)+'_'+str(p1)+'.pdf',rasterized=True)
plt.savefig('fisher_barnard_sim'+str(fixn1)+'_'+str(n)+'_'+str(pn1)+'_'+str(p1)+'.png')
plt.show()
