import math
import numpy as np
import scipy.stats as stats
import copy

def binomial(n , k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def logprob(x , a , pi):
    x = float(x)
    a = float(a)
    pi = float(pi)
    return binomial(n - n1 , x - a) + binomial(n1 , a) + x * math.log(pi) + (n - x) * math.log(1.0 - pi)

def sumlogs(a , b):
    return max(a,b) + np.log1p(math.exp(min(a,b)-max(a,b)))

def pre_load_pvalue_barnard(x , pi):
    p_values = list()
    p_values_log = list()
    for x_ in range(1,int(n)):
        alb = max(0,n1-(n-x_))
        aub = min(n1,x_)
        for a_ in range(int(alb),int(aub)+1):
            p = logprob(x_ , a_ , pi)
            p_values.append(math.exp(p))
            p_values_log.append(p)
    p_values.sort()
    p_values_log.sort()
    #print(p_values_log)
    return np.array(p_values) , np.array(p_values_log)

def pvalue_barnard_preloaded(x , a , pi , p_values , p_values_log):
    pval = 0
    p0 = math.exp(logprob(x , a , pi))
    pval = np.sum( p_values[ p_values <= p0] )
    return pval

def pvalue_barnard(a0 , a1 , n1 , n0 , n , pi):
    pval = -99999
    x = a0 + a1
    a = a1
    p0 = logprob(x , a , pi)
    for x_ in range(1,int(n)):
        alb = max(0,n1-(n-x_))
        aub = min(n1,x_)
        for a_ in range(int(alb),int(aub)+1):
            p = logprob(x_ , a_ , pi)
            if p <= p0:
                pval = sumlogs(pval , p)
    return pval


def pvalue_fisher(x , a):
    a0 = x - a
    oddsratio, pvalue_f = stats.fisher_exact([[a, n1-a], [a0, n0-a0]])
    return pvalue_f

n = 10000.
pn1 = 0.25
n1 = n* pn1
n0 = n - n1
p1 = 0.1
x = p1 * n
max_a = min(n1 , x)
min_a = max(0 , x - (n - n1))
pi = x / n

x_ = list()
y_ = list()
y_f = list()
y_b = list()

p_values, p_values_log = pre_load_pvalue_barnard(x , pi)
print("done preloading")
p_values_barnard_opt = list()
as_computed = list()

mid_a = int(float(x) * pn1)

a = mid_a
pvalue_f = 1.0
min_pvalue_ = 1e-10
while pvalue_f > min_pvalue_ and a <= max_a:
    print(a)
    pvalue_f = pvalue_fisher(int(x) , a)
    if pvalue_f > min_pvalue_:
        pvalue_b = pvalue_barnard_preloaded(x , float(a) , pi , p_values, p_values_log)
        #pvalue_b = math.exp(pvalue_barnard(x , float(a) , pi))
        y_.append(pvalue_f / pvalue_b)
        y_f.append(pvalue_f)
        y_b.append(pvalue_b)
        x_.append(a)
        out_ = str(a) + " " + str(pvalue_f)+ " " + str(pvalue_b)+ " " + str(pvalue_f / pvalue_b)
    else:
        out_ = str(a) + " " + str(pvalue_f)
    a = a + 1

a = mid_a-1
pvalue_f = 1.0
while pvalue_f > min_pvalue_ and a >= min_a:
    print(a)
    pvalue_f = pvalue_fisher(int(x) , a)
    if pvalue_f > min_pvalue_:
        pvalue_b = pvalue_barnard_preloaded(x , float(a) , pi , p_values, p_values_log)
        #pvalue_b = math.exp(pvalue_barnard(x , float(a) , pi))
        y_.append(pvalue_f / pvalue_b)
        y_f.append(pvalue_f)
        y_b.append(pvalue_b)
        x_.append(a)
        out_ = str(a) + " " + str(pvalue_f)+ " " + str(pvalue_b)+ " " + str(pvalue_f / pvalue_b)
    else:
        out_ = str(a) + " " + str(pvalue_f)
    a = a - 1

p_values_barnard_opt = list()

for index, a in enumerate(x_):
    p_values_barnard_opt.append(0.0)




pis = np.linspace(pi*0.9 , pi*1.1 , 100)
for jindex , pi in enumerate(pis):
    p_values, p_values_log = pre_load_pvalue_barnard(x , pi)
    print("done preloading "+str(jindex+1)+"/"+str(len(pis)))
    for index, a in enumerate(x_):
        pvalue_b_o = pvalue_barnard_preloaded(x , float(a) , pi , p_values, p_values_log)
        p_values_barnard_opt[index] = max(pvalue_b_o , p_values_barnard_opt[index])
        out_ = str(a) + " " + str(p_values_barnard_opt[index])

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
f, ax1 = plt.subplots(figsize=(5, 4))
ax1.set_yscale('log', nonposy="clip")
ax1.set_xscale('log', nonposx="clip")

maxy_ = float(np.max(y_))

pi = x / n
max_value = x_[np.argmax(y_)]
print(np.min(y_))
#ax1.axvline(x=pi , color="black" , linewidth = 1)
#ax1.axvline(x=max_value , color="red" , linewidth = 1)
#ax1.axhline(y=0.05 , color="red" , linewidth = 1)
#ax1.axhline(y=0.1 , color="red" , linewidth = 1)
#ax1.set_xlim(x - dev, x + dev)
#ax1.set_ylim(10**-5, 1.0)
plt.ylabel(r'Test $p$-value')
plt.xlabel(r'Exact $p$-value')
#plt.title(r'$x = $'+str(int(x))+' $n = $'+str(int(n))+' $n_1 = $'+str(int(n1)))
plt.plot(y_b , y_f , '.', c="red" , markersize=2 , label=r'Fisher')
plt.plot(y_b , p_values_barnard_opt , '.', c="blue" , markersize=2 , label=r'Barnard')
min_ = min( min(y_b , y_f) ) *0.95
max_ = max( max(y_b , y_f) ) *1.05
plt.plot([min_ , max_] , [min_ , max_] , linewidth=0.5 , c="black")
ax1.set_xticks(ax1.get_yticks())
ax1.set_yticks(ax1.get_yticks())
labels_y = np.linspace( min( min(y_b , y_f) ) , max_y , 5)
#plt.yticks( labels_y , labels_y )
#ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3e'))
ax1.set_ylim(max_, min_)
ax1.set_xlim(max_, min_)
plt.legend()
minor_ticks = [10**-2 , 10**-4, 10**-6, 10**-8]
print(minor_ticks)
ax1.set_xticks(minor_ticks, minor=True)
ax1.set_yticks(minor_ticks, minor=True)
plt.grid(b=True, which='both', color='0.65', linestyle='-')
f.patch.set_facecolor('xkcd:white')
plt.tight_layout()
plt.savefig("pvalues_fb_comp.pdf")
plt.show()
