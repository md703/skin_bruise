#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 11:10:15 2020

@author: md703
"""

import numpy as np
import matplotlib.pyplot as plt

def export_musp(wl, A, K):
    musp = 1000 * A * wl**(-K)
    return musp

pl_folder = "johnson_pathlength"
wl = np.arange(450, 600+1, 5)[:, None]
A = np.array([1980000, 350322])
K = np.array([2.833, 2.631])

musp = export_musp(wl, A, K)

plt.plot(wl, 1.3*musp[:, 0], label="stratum corneum")
plt.plot(wl, musp[:, 0], label="epidermis")
plt.plot(wl, musp[:, 1], label="dermis")
plt.legend()
plt.xlabel("wavelength [nm]")
plt.ylabel("musp [1/cm]")
plt.title("musp of each skin layer")
plt.savefig("{}/musp of each skin layer".format(pl_folder), dpi=300)