# python modules
import pandas as pd
import numpy as np
import re
from collections import defaultdict
from scipy.signal import medfilt
from scipy.signal import argrelextrema
from math import isclose
import glob

import os
import seaborn as sns
import matplotlib.pyplot as plt



def get_local_global_alpha_value(F, D_local_global):
    all_dfs = []

    for f in F:
        df = pd.read_csv(f, sep="\t", index_col=0)    

        if df.shape[0] == 0:
            continue

        else:
            all_dfs.append(df)
    
    
    if len(all_dfs) == 0:        
        
        a = {'2nd-eigenvalue': np.nan,
             '2nd-largest-cc-size': np.nan,
             'almost-disconnected-component-count': np.nan,
             'alpha': np.nan,
             'connected-component-count': np.nan,
             'density': np.nan,
             'edge-count': np.nan,
             'largest-cc-size': np.nan,
             'mean-k': np.nan,
             'vertex-count': np.nan}
        
        D_local_global["alpha_max"] = a
        D_local_global["alpha_exist"] = a
        return pd.DataFrame()
    
    
    
    df = pd.concat(all_dfs).groupby("alpha").max() #skipna=True)    
    df.reset_index(inplace=True)

    max_ac = df[(df["2nd-eigenvalue"] < 1) & \
                             (df["almost-disconnected-component-count"] > 1)]
    if max_ac.shape[0] > 0:
        row_alpha_exist = max_ac[max_ac["alpha"] == max_ac["alpha"].min()]
        
        row_alpha_max_ = max_ac[max_ac["almost-disconnected-component-count"] == max_ac["almost-disconnected-component-count"].max()]
        row_alpha_max  = row_alpha_max_[row_alpha_max_["alpha"] == row_alpha_max_["alpha"].min()]
        
        D_local_global["alpha_max"] = row_alpha_max.squeeze().to_dict()
        D_local_global["alpha_exist"] = row_alpha_exist.squeeze().to_dict()
    else:
        a = {'2nd-eigenvalue': np.nan,
             '2nd-largest-cc-size': np.nan,
             'almost-disconnected-component-count': np.nan,
             'alpha': np.nan,
             'connected-component-count': np.nan,
             'density': np.nan,
             'edge-count': np.nan,
             'largest-cc-size': np.nan,
             'mean-k': np.nan,
             'vertex-count': np.nan}
        
        D_local_global["alpha_max"] = a
        D_local_global["alpha_exist"] = a
        
    return df

def get_significance_t_values(F, D, alpha=0.05, min_power=0.8):
    
    all_power_df = []
    for f in F:
        with open(f, "r") as f_in:
            line1 = f_in.readline()
            alpha, sample_size, r = re.findall("^\# Critical Pearson correlation for obtaining alpha significance of (\d*\.?\d*) with sample size of (\d*) is (\d*\.?\d*)$", line1)[0]
            D["TypeI-"+ str(alpha)] = float(r)

        df = pd.read_csv(f, sep="\t", skiprows=2, index_col=0)
        if df.shape[0] == 0:
            continue
        
        else:
            all_power_df.append(df)
    
    if len(all_power_df) == 0:
        D["TypeI-"+ str(alpha)] = np.nan      
        D["Power-"+ str(min_power)] = np.nan
        return pd.DataFrame()
    
    power_df = pd.concat(all_power_df).groupby("r").max()#skipna=True)    
    power_df.reset_index(inplace=True)
    
    D["Power-"+ str(min_power)] = power_df[power_df["power"] >= min_power]["r"].min()
    
    return power_df

def get_interative_t_values(F, D, d_min_t={"general":0}):
    
    all_dfs = []
    
    for f in F:
        df = pd.read_csv(f, sep="\t")    
        
        if df.shape[0] == 0:
            continue
        
        else:
            all_dfs.append(df)
        
    if len(all_dfs) == 0:
        return pd.DataFrame()
    
    
    df = pd.concat(all_dfs).groupby("threshold").max()# skipna=True)    
    df.sort_values("threshold", inplace=True)
    df.reset_index(inplace=True)
    df = df[df["threshold"] > d_min_t["general"]]
    df.reset_index(inplace=True)

    D['aplestin'] = np.nan
    D['cc_inflection'] = np.nan
    D['density_min'] = np.nan
    D['elo_clustering'] = np.nan
    D['gupta_clustering'] = np.nan
    D['mcr_2'] = np.nan
    D['mcr_3'] = np.nan
    D['mcr_max'] = np.nan
    D['rmt'] = np.nan
    D['scale_free'] = np.nan
    D['single_component'] = np.nan
    D['spectral_methods'] = np.nan
    D['whole_graph'] = np.nan        
    D['Power-0.8'] = np.nan        
    D['TypeI-0.01'] = np.nan        

    if df.shape[0] <3:
        return df
    
    # single component
    # First point before more than one cc appears
    D['single_component'] = df[df["connected-component-count"] == 1]["threshold"].max()
    
    # whole graph
    D['whole_graph'] = df[df["vertex-count"] ==  df["vertex-count"].max()]["threshold"].max()

    # giant cc down inflection
    diffs = np.diff(df["largest-cc-size"])
    drop_cutoff = diffs.std()
    diffs_b = (abs(diffs) > drop_cutoff) & (diffs < 0)
    if sum(diffs_b) > 0:
        D['cc_inflection'] = df["threshold"][np.argmax(diffs_b) + df.index[0]]
        
    # density minimum
    D['density_min'] = df[df["density"] ==  df["density"].min()]["threshold"].min()
    
    # scale free - lowest value where KS-pvalue >0.1
    D['scale_free'] = df[df["scale-free-KS-p-value"] > 0.1]["threshold"].min()

    # spectral methods lowest t with largest number 
    # of almost-disconnected-components
    # check if there are
    if "spectral" in d_min_t.keys():
        subdf = df[df["threshold"] > d_min_t["spectral"]]
        subdf.reset_index(inplace=True)
    else: 
        subdf = df
                   
    subdf.loc[subdf["almost-disconnected-component-count"] == -1, \
           "almost-disconnected-component-count"] = np.nan
    if not np.all(np.isnan(subdf["almost-disconnected-component-count"])):
        # posiblilites: fileder < 2 (?)
        max_ac = subdf[subdf["2nd-eigenvalue"] < 2]
        if max_ac.shape[0] > 0:
            D['spectral_methods'] = max_ac[max_ac["almost-disconnected-component-count"] == max_ac["almost-disconnected-component-count"].max()]["threshold"].min()

    # RMT somewhere between 
    # poisson chi pvalue is > 0.05 (i.e. Poisson), 
    # goe chi2 pvalue > 0.05 (i.e. not GOE)
    if "rmt" in d_min_t.keys():
        subdf = df[df["threshold"] > d_min_t["rmt"]]
        subdf.reset_index(inplace=True)
    else: 
        subdf = df                   
                   
    D['rmt'] = (subdf[subdf["poisson-pvalue"] > 0.05]["threshold"].min() + \
                 subdf[subdf["goe-pvalue"] < 0.05]["threshold"].max()) /2

    # Maximal clique ratio
    if "mcr" in d_min_t.keys():
        subdf = df[df["threshold"] > d_min_t["mcr"]]
        subdf.reset_index(inplace=True)
    else: 
        subdf = df                      
                   
                   
    subdf.loc[subdf["maximal-clique-count"] == -1, "maximal-clique-count"] = np.nan
    if not np.all(np.isnan(subdf["maximal-clique-count"])):
        subdf["maximal-clique-ratio"] = subdf["maximal-clique-count"]/np.append(subdf["maximal-clique-count"][1:].values, (np.nan))
        subdf['maximal-clique-ratio'].replace(np.inf, np.nan, inplace=True)
        D['mcr_2'] = subdf[subdf["maximal-clique-ratio"] >= 2]["threshold"].min()
        D['mcr_3'] = subdf[subdf["maximal-clique-ratio"] >= 3]["threshold"].min()
        D['mcr_max'] =  subdf[subdf["maximal-clique-ratio"] == subdf["maximal-clique-ratio"].max()]["threshold"].min()

    # apelstin
    Nsv = df["edge-count"]/df["vertex-count"]
    dNsv_dt = np.gradient(Nsv, df["threshold"])
    df["dNsv_dt"] = dNsv_dt
    D['aplestin'] = df["threshold"][dNsv_dt > 0].min()

    if not np.all(np.isnan(df["clustering-coefficient"])):
        # gupta clustering coefficient
        # largest t with sharp increase
        # estimate by first point where at least 
        # 3 differences are larger than 0.5 * stddev
        found_gupta = False
        diffs = np.diff(df["clustering-coefficient"])
        if (np.sum(diffs>0) > 3): # first check: is there consistent increase in clustering coeficeint?
            cutoff = np.nanstd(diffs[diffs>0])*0.3
            prev_prev_d = diffs[0]
            prev_d = diffs[1]
            i = 0
            for d in diffs[2:]:
                if np.all(np.array([prev_prev_d, prev_d, d]) > cutoff) and (i > 1):
                    found_gupta = True
                    break
                prev_prev_d = prev_d
                prev_d = d
                i += 1

        if found_gupta:
            D['gupta_clustering'] = df["threshold"][i]
        
        if not np.all(np.isnan(df["random-clustering-coefficient"])):
            # elo clustering coefficient
            # first local maximum
            C0_diffs =  medfilt(df["clustering-coefficient"] - df["random-clustering-coefficient"])
            df["elo-clustering-coefficient-diffs"] = C0_diffs
            D['elo_clustering'] = df["threshold"][argrelextrema(C0_diffs, np.greater_equal)[0][0] + df.index[0]].min()

    return df


def get_result_per_prefix(prefix):
    iterative_results = glob.glob(prefix + "*.iterative.txt")
    significance_results = glob.glob(prefix + "*.statistical_errors.txt")
    local_global_results = glob.glob(prefix + "*.local_global.txt")
    
    D = {}
    alpha = np.nan
    for method, F in [("iterative_result", iterative_results), 
                  ("significance_result", significance_results), 
                  ("local_global_result", local_global_results)]:
        if len(F)>0:
            if method == "iterative_result":
                df = get_interative_t_values(F, D)

            elif method == "significance_result":
                power_df = get_significance_t_values(F, D, min_power=0.8)

            elif method == "local_global_result":
                df, alpha = get_local_global_alpha_value(F)
    return D, alpha


def plot_t_vs_EV(df_plot, D):
    
    annotations = {}
    labels = defaultdict(set)

    for k1, t1 in D.items():
        if not np.isnan(t1):
            is_in = False
            for k2, t2 in annotations.items():
                if isclose(t1, t2):
                    labels[k2].update([k1])
                    annotations[k2] = np.mean([t1, t2])
                    is_in = True
            if not is_in:
                labels[k1].update([k1])
                annotations[k1] = t1

    
    # number vertices and number edges vs thresholds
    with sns.plotting_context("paper"):
        fig, ax = plt.subplots(figsize=(12, 6)) # long, high

        ax.plot(df_plot["vertex-count"], alpha=0.6, label="Vertex count")

        ax.plot(np.nan, np.nan, alpha=0.6, label="Edge-count", color="orange")

        ax_twin = ax.twinx()
        ax_twin.plot(df_plot["edge-count"], alpha=0.6, color="orange", label="Edge count")

        xmin, xmax, ymin, ymax = ax.axis()
        y_txt = (ymax - ymin)/100 + 10

        for key, value in annotations.items():
            ax.plot(value, y_txt, 
                    marker=11, #"CARETDOWNBASE"
                    markersize=10,
                    color="red")

            txt = " &\n".join(labels[key])
            ax.annotate(txt, (value, y_txt), 
                        rotation=-45, 
                        horizontalalignment='right', 
                        verticalalignment='bottom')

        ax.set(xlabel="Threshold")
        ax.set(ylabel='Vertex count')
        ax_twin.set(ylabel="Edge count")
        ax_twin.grid(False)

        ax.legend()
        ax.legend(loc=1)
        sns.despine(ax=ax,      right=True, left=False, bottom=True, top=True)

        plt.tight_layout()
        
