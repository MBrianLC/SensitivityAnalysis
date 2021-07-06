import numpy as np
import pandas as pd
import scipy
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

namesSize = ['1MB','2MB','3MB','4MB','5MB','6MB','7MB','8MB','9MB','10MB']
namesBias = ['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1']

# It is used to color the sensitivity matrix. It uses three levels of significance (0.1, 0.5, and 0.01)
def frontier_color(col):
    df = col.copy()
    df.loc[:] = ''
    v = [0.001,0.01,0.05]
    for i in range(10):
        if col[i] < v[0] or 1 - col[i] < v[0]:
            df.iloc[i] = 'color: red'
        elif col[i] < v[1] or 1 - col[i] < v[1]:
            df.iloc[i] = 'color: orange'
        elif col[i] < v[2] or 1 - col[i] < v[2]:
            df.iloc[i] = 'color: #cccc00'
        else:
            df.iloc[i] = 'color: green'
    return df

# Same as the previous method, but with only one level of significance (0.01)
def frontier_color_2(col):
    df = col.copy()
    df.loc[:] = ''
    v = 0.01
    for i in range(10):
        if col[i] < v or 1 - col[i] < v:
            df.iloc[i] = 'color: red'
        else:
            df.iloc[i] = 'color: green'
    return df

# Applies the coloring methods to the sensitivity matrix
def frontier(df):
    return df.style.apply(frontier_color,axis=0)

# Kolmogorov-Smirnov test on the list p[name]    
def ksTest(p,name):
    return stats.kstest(p[name],'uniform')[0]

# Applies the Kolmogorov-Smirnov test on blocks of 10 elements from the list p[name]
def ksTestAll(p,name):
    ks = [[] for i in range(10)]
    for i in range(10):
        for j in range(10):
            ks[i].append(ksTest(p[10*(10*i+j):10*(10*i+j+1)],name))
    dks = pd.DataFrame(ks)
    dks.columns = namesBias
    dks.index = namesSize
    return frontier(dks)

# Proportion of failed tests of the p-values list p
def ftests(p,alpha):
    return sum(p < alpha)/p.size

# Proportion of failed tests on blocks of 10 elements from the list p[name]
def ftAll(p,name):
    ft = [[] for i in range(10)]
    for i in range(10):
        for j in range(10):
            ft[i].append(ftests(np.array(p[10*(10*i+j):10*(10*i+j+1)][name]).transpose(),0.01))
    dft = pd.DataFrame(ft)
    dft.columns = namesBias
    dft.index = namesSize
    return frontier(dft)

# Representation of the sensitivity frontier    
def pltStep(v1,v2,v3,name):
    x = np.linspace(0,10,11)
    y1 = np.array(v1)
    y2 = np.array(v2)
    y3 = np.array(v3)

    fig, ax = plt.subplots()
    ax.plot(x, y1, 'b', linewidth=2)
    ax.plot(x, y2, 'g', linewidth=2)
    ax.plot(x, y3, 'r', linewidth=2)
    ax.set_ylim(bottom=0)

    ix = np.linspace(0,10,11)
    iy1 = y1
    verts = [(0, 0), *zip(ix, iy1), (10, 0)]
    poly1 = Polygon(verts, facecolor='0.9', edgecolor='0.2')
    ax.add_patch(poly1)

    iy2 = y2
    verts = [(0, 0), *zip(ix, iy2), (10, 0)]
    poly2 = Polygon(verts, facecolor='0.75', edgecolor='0.2')
    ax.add_patch(poly2)

    iy3 = y3
    verts = [(0, 0), *zip(ix, iy3), (10, 0)]
    poly3 = Polygon(verts, facecolor='0.6', edgecolor='0.2')
    ax.add_patch(poly3)

    ax.xaxis.set_ticks_position('bottom')

    ax.set_xticks(x)
    ax.set_xticklabels(['0']+namesBias)
    ax.set_yticks(x[1:])
    ax.set_yticklabels(namesSize[::-1])
    plt.title(name)
    plt.show()

# Numerical sensitivity of v1, v2 and v3 (values of the sensitivity frontier for the three biased generators). f is the numerical integration method 
def sensf(f,v1,v2,v3,testNames,nTests)    
    nTests = 5
    k = [[] for i in range(nTests)]
    for i in range(nTests):
        k[i].append(f(v1[i]))
        k[i].append(f(v2[i]))
        k[i].append(f(v3[i]))
    pk = pd.DataFrame(np.array(k))
    pk.index = testNames
    pk.columns = ["ε-Hole","σ-Counter","T-Counter"]
    return pk

# Numerical sensitivity of the vectors v1, v2 and v3 (using scipy.integrate.trapz as the numerical integration method)
def sensitivity(v1,v2,v3,testNames,nTests)
    return sensf(scipy.integrate.trapz,v1,v2,v3,testNames,nTests)

# Percentage of failed tests of the set of lists of p-values p (where alpha is the level of significance)
def failedtests(p,alpha,nTests):
    nTests = len(p)
    ft = np.zeros(nTests)
    for i in range(nTests):
        ft[i] = sum(p[i] < alpha)/p[i].size
    return ft

# Percentage of failed tests of the three sets of lists of p-values p (one for each biased generator)
def failedTests(pE,pS,pT,testNames,nTests)
    ft = []
    ft.append(failedtests(np.array(pE).transpose(),0.01,nTests))
    ft.append(failedtests(np.array(pS).transpose(),0.01,nTests))
    ft.append(failedtests(np.array(pT).transpose(),0.01,nTests))
    dft = pd.DataFrame(ft)
    dft.index = ["ε-Hole","σ-Counter","T-Counter"]
    dft.columns = testNames
    return dft