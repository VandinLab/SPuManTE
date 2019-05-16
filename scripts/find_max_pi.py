import math
import numpy as np

def binomial(n , k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def logprob(x , a , pi):
    x = float(x)
    a = float(a)
    pi = float(pi)
    return binomial(n - n1 , x - a) + binomial(n1 , a) + x * math.log(pi) + (n - x) * math.log(1.0 - pi)

def sumlogs(a , b):
    return max(a,b) + np.log1p(math.exp(min(a,b)-max(a,b)))

def pvalue_barnard(x , a , pi):
    pval = -99999
    p0 = logprob(x , a , pi)
    for x_ in range(1,int(n)):
        alb = max(0,n1-(n-x_))
        aub = min(n1,x_)
        for a_ in range(int(alb),int(aub)+1):
            p = logprob(x_ , a_ , pi)
            if p <= p0:
                pval = sumlogs(pval , p)
    return pval



n = 1000.0
n1  = 500.0
x = 100.0
a = 40.0
num_points = 1000
grains = [0.1]

for grain in grains:

    y_ = list()
    x_ = list()

    pi = x / n
    pis = np.linspace(pi*(1.0-grain),pi*(1.0+grain),num_points)

    pis = list(pis)
    pis.append(pi)
    pis.sort()
    #print pis

    pi = x / n
    simplepval = pvalue_barnard(x , a , pi)
    print simplepval
    print math.exp(simplepval)

    ids = list()
    i = 0
    from multiprocessing import Pool
    pool = Pool()
    for pi in pis:
        ids.append( pool.apply_async(pvalue_barnard , [x , a , pi]))
        #print "p"+str(i)
        i = i+1
    i = 0
    for pi in pis:
        x_.append(pi)
        y_.append(math.exp(ids[i].get()))
        print "*p"+str(i)
        i = i+1

    # for pi in pis:
    #     x_.append(pi)
    #     y_.append(math.exp(pvalue_barnard(x , a , pi)))

    #print y_



    import matplotlib as mpl
    mpl.use('Agg')
    mpl.rcParams['figure.dpi'] = 300
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mtick
    f, ax1 = plt.subplots(figsize=(6, 3))
    ax1.set_yscale('log', nonposy="clip")

    pi = x / n
    max_value = x_[np.argmax(y_)]
    maxy_ = float(np.max(y_))
    ax1.axhline(y=math.exp(simplepval) , color="black" , linewidth = 1)
    ax1.axhline(y=maxy_ , color="red" , linewidth = 1)
    ax1.axvline(x=pi , color="black" , linewidth = 1)
    ax1.axvline(x=max_value , color="red" , linewidth = 1)
    #ax1.set_xlim(x - dev, x + dev)
    #ax1.set_ylim(10**-5, 1.0)
    plt.ylabel(r'$p$-value')
    plt.xlabel(r'$\pi$')
    plt.plot(x_ , y_ , '.', markersize=6)
    min_y = float(min(y_))*0.95
    max_y = float(max(y_))*1.05
    print str(min_y)
    labels_y = np.linspace( min_y , max_y , 5)
    plt.yticks( labels_y , labels_y )
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.3e'))
    ax1.set_ylim(min_y, max_y)
    ax1.set_xlim(pi*(1. - (1.*grain + grain/num_points)), pi*(1. + (1.*grain + grain/num_points)))
    plt.tight_layout()
    plt.savefig("pis_values_"+str(grain)+".pdf")
    plt.show()

    print pi
    print x_[np.argmax(y_)]
    print x_[np.argmax(y_)]*n
    print np.max(y_)
    print math.exp(simplepval)
    print np.max(y_) / math.exp(simplepval)
