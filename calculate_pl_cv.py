#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 16:18:04 2020

@author: md703
"""

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

#%% paramters
detector_num = 3
layer_num = 3
run_num = 10
wl_file = "sim_wl.txt"
pl_folder = "johnson_pathlength"
pl_type = "johnson_mbll"
pl_file = "average_pathlength.txt"
reflectance_file = "GPUMC_output.txt"

#%% main
# read wl values
wl = pd.read_csv("{}/{}_0/{}".format(pl_folder, pl_type, wl_file), 
                 sep='\t   ', engine='python')
wl = wl.columns.values.astype(float).astype(int)
wl_num = wl.size

# create pl, reflectance container
all_wl_reflectance = {}
all_wl_pathlength = {}

for detector_idx in range(1, detector_num+1):
    all_wl_reflectance["detector_{}".format(detector_idx)] = {}
    all_wl_pathlength["detector_{}".format(detector_idx)] = {}
    for wl_idx in range(1, wl_num+1):
        all_wl_reflectance["detector_{}".format(detector_idx)]["wl_{}".format(wl_idx)] = {}
        all_wl_pathlength["detector_{}".format(detector_idx)]["wl_{}".format(wl_idx)] = {}

# read pl, reflectance values and calculate ...
reflectance_mean_set = np.empty(0)
reflectance_cv_distribution = np.empty(0)
pl_mean_set = np.empty(0)
pl_cv_distribution = np.empty(0)
monitor_info = "detector_{}, {}, {}nm ------> reflectance cv: {:.3%}, mean: {:.4f} ; pathlength cv: {:.3f}, mean: {:.4f}"

for detector_idx in range(1, detector_num+1):
    
    for wl_idx in range(1, wl_num+1):
        wl_dir = "wl_{}".format(wl_idx)
        reflectance_set = np.empty(run_num)
        pl_set = np.empty(run_num)
        
        for run_idx in range(run_num):
            run_dir = "{}_{}".format(pl_type, run_idx)
            
            # read reflectance
            reflectance_temp = pd.read_csv(os.path.join(pl_folder, run_dir, wl_dir, reflectance_file), 
                                           sep='\s',
                                           engine='python')
            reflectance_temp = reflectance_temp.columns.values.astype(np.float64)
            # get reflectance at each sds(detector)
            reflectance_set[run_idx] = reflectance_temp[detector_idx-1]
            
            # read pl
            pl_temp = pd.read_csv(os.path.join(pl_folder, run_dir, wl_dir, pl_file), sep='\t')
            pl_temp = pl_temp.drop(columns='Unnamed: 9', axis=1).columns.values.astype(np.float64)
            pl_temp = pl_temp.reshape(detector_num, layer_num)            
            # get pathlength in dermis at each sds(detector)
            pl_set[run_idx] = pl_temp[detector_idx-1][-1]
            
        # get reflectance statistics-related info
        all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["values"] = reflectance_set
        all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["mean"] = reflectance_set.mean()
        all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["std"] = reflectance_set.std(ddof=1)
        all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["cv"] = reflectance_set.std(ddof=1) / reflectance_set.mean()
        
        # get pl statistics-related info
        all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["values"] = pl_set
        all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["mean"] = pl_set.mean()
        all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["std"] = pl_set.std(ddof=1)
        all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["cv"] = pl_set.std(ddof=1) / pl_set.mean()
        
        print(monitor_info.format(detector_idx, wl_dir, wl[wl_idx-1],
                                  all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["cv"],
                                  all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["mean"],
                                  all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["cv"],
                                  all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["mean"]))
        
        reflectance_mean_set = np.append(reflectance_mean_set, all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["mean"])
        reflectance_cv_distribution = np.append(reflectance_cv_distribution, all_wl_reflectance["detector_{}".format(detector_idx)][wl_dir]["cv"])
        
        pl_mean_set = np.append(pl_mean_set, all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["mean"])
        pl_cv_distribution = np.append(pl_cv_distribution, all_wl_pathlength["detector_{}".format(detector_idx)][wl_dir]["cv"])
        
#%% plot mean of 10 reflectance of each wl
reflectance_mean_set = reflectance_mean_set.reshape(detector_num, wl_num)
reflectance_dict = {"wl":wl}
for detector_idx in range(1, detector_num+1):
    reflectance_dict["detector_{}_reflectance_mean".format(detector_idx)] = reflectance_mean_set[detector_idx-1]
reflectance_df = pd.DataFrame(reflectance_dict, columns=['wl'] + ["detector_{}_reflectance_mean".format(detector_idx) for detector_idx in range(1, detector_num+1)])
# np.savetxt(r'{}/pathlength.txt'.format(pl_folder), pl_df.values)
# pl_df.to_csv('{}/pathlength.csv'.format(pl_folder))

#%% save mean of 10 pathlength of each wl
pl_mean_set = pl_mean_set.reshape(detector_num, wl_num)
pl_dict = {"wl":wl}
for detector_idx in range(1, detector_num+1):
    pl_dict["detector_{}_pl_mean".format(detector_idx)] = pl_mean_set[detector_idx-1]
pl_df = pd.DataFrame(pl_dict, columns=['wl'] + ["detector_{}_pl_mean".format(detector_idx) for detector_idx in range(1, detector_num+1)])
# np.savetxt(r'{}/pathlength.txt'.format(pl_folder), pl_df.values)
# pl_df.to_csv('{}/pathlength.csv'.format(pl_folder))

#%% plot reflectance cv distribution
reflectance_cv_distribution = reflectance_cv_distribution.reshape(detector_num, wl_num)
plt.figure(figsize=(10, 5))

for detector_idx in range(1, detector_num+1):
    plt.plot(wl, reflectance_cv_distribution[detector_idx-1], marker='.', label="detector_{}".format(detector_idx))
    
plt.axhline(y=0.01, color='r')
plt.legend()
plt.xlabel("wavelength [nm]")
plt.ylabel("CV [-]")
plt.title("reflectance CV distribution for {} detector(s)".format(detector_num))
# plt.savefig("{}/pathlength_cv_distribution".format(pl_folder), dpi=300, bbox_inches='tight')
plt.show()

#%% plot reflectance distribution
for detector_idx in range(1, detector_num+1):
    plt.plot(reflectance_df["wl"].values, reflectance_df["detector_{}_reflectance_mean".format(detector_idx)].values, marker='.', 
             label="detector_{}".format(detector_idx))
plt.legend()
plt.xlabel("wavelength [nm]")
plt.ylabel("reflectance [-]")
plt.title("reflectance_distribution for {} detector(s)".format(detector_num))
# plt.savefig("{}/pathlength_distribution".format(pl_folder), dpi=300, bbox_inches='tight')
plt.show()

#%% plot reflectance cv distribution
pl_cv_distribution = pl_cv_distribution.reshape(detector_num, wl_num)
plt.figure(figsize=(10, 5))

for detector_idx in range(1, detector_num+1):
    plt.plot(wl, pl_cv_distribution[detector_idx-1], marker='.', label="detector_{}".format(detector_idx))
    
plt.axhline(y=0.01, color='r')
plt.legend()
plt.xlabel("wavelength [nm]")
plt.ylabel("CV [-]")
plt.title("pathlength CV distribution for {} detector(s)".format(detector_num))
# plt.savefig("{}/pathlength_cv_distribution".format(pl_folder), dpi=300, bbox_inches='tight')
plt.show()

#%% plot pathlength distribution
for detector_idx in range(1, detector_num+1):
    plt.plot(pl_df["wl"].values, pl_df["detector_{}_pl_mean".format(detector_idx)].values, marker='.', 
             label="detector_{}".format(detector_idx))
plt.legend()
plt.xlabel("wavelength [nm]")
plt.ylabel("pathlength [cm]")
plt.title("pathlength_distribution for {} detector(s)".format(detector_num))
# plt.savefig("{}/pathlength_distribution".format(pl_folder), dpi=300, bbox_inches='tight')
plt.show()




