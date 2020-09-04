#! /usr/bin/env python3
# -*- Coding: UTF-8 -*-

r"""
Coefficient of Distribution in Analysical Chemistry of weak acid.
"""

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Ka = [None, ]
# Read the Ka from a excel 'Ka.xlsx', or from command arguments. You can only change the code here.
if len(argv) > 1:
    Ka.extend(list(map(float, argv[1:])))
else:
    Ka.extend(list(pd.read_excel('Ka.xlsx', header = None, squeeze = True))) 
ndiv = 10 # steps when pH changed 1

ngrp = len(Ka)
nhyd = ngrp - 1
Ign1 = lambda _: str(_) if _ != 1 else str()
npoints = ndiv * 14 + 1
x = np.linspace(0, 14, npoints)
H = np.power(10, -x)
y = np.ones((ngrp, npoints), dtype = np.float64)
for i in np.arange(ngrp):
    for k in 1 + np.arange(i): y[i] *= Ka[k]
    y[i] *= H ** (nhyd - i)
s = sum(y)
y /= s

fig, ax = plt.subplots()
for i in range(ngrp):
    if i == 0:
        lab = '$\\mathrm{H_{%s}A}$' % Ign1(nhyd)
    elif i == nhyd:
        lab = '$\\mathrm{A^{%s-}}$' % Ign1(nhyd)
    else:
        lab = '$\\mathrm{H_{%s}A^{%s-}}$' % (Ign1(nhyd - i), Ign1(i))
    ax.plot(x, y[i], label = lab)
ax.legend()
ax.set_xticks(np.linspace(0, 14, 15))
ax.set_yticks(np.linspace(0, 1, 11))
ax.set_xlabel('pH')
ax.set_ylabel('Distribution coefficient')
ax.set_title('Distribution Coefficient with pH')
fig.savefig('Distribution_Coefficient_with_pH.png')
plt.show()

