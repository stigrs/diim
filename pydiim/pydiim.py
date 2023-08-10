# Copyright (c) 2023 Stig Rune Sellevag
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

"""Provides a Python wrapper for the DIIM C++ code."""

import numpy as np
import pandas as pd
import json
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import AutoLocator, AutoMinorLocator


def score_to_interdependency(score_file, amat_file, score_scale):
    """Transform score values to interdependencies."""
    try:
        subprocess.run(["diim_gen", score_file, amat_file, str(score_scale)], check=True)
    except subprocess.CalledProcessError as exc:
        print(f"{exc}")
    except Exception as exc:
        print(f"{exc}")


def read_csv(filename, encoding="latin1"):
    """Helper function for reading CSV files."""
    return pd.read_csv(filename, encoding=encoding)


def excel_writer(filename):
    """Helper function for writing to Excel XLSX file."""
    return pd.ExcelWriter(filename, engine="xlsxwriter")


def bar_plot(xdata, ydata, xlabel=None, ylabel=None, title=None, figsize=(10, 5), dpi=300):
    """Helper function for creating IIM plots."""
    _, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.bar(xdata, ydata)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    plt.xticks(rotation=90)


def grouped_bar_plot(xtick_labels, 
                     data, 
                     xlabel=None, 
                     ylabel=None, 
                     legend=None, 
                     title=None, 
                     bar_width=0.4, 
                     figsize=(10, 5), 
                     dpi=300):
    """Helper function for creating grouped IIM bar plots."""
    _, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.yaxis.grid(color='gray', linestyle='dashed')

    x_pos = np.arange(len(xtick_labels))

    # Loop over data:
    for i, (group, values) in enumerate(data.items()):
        pos = x_pos + (i * bar_width)
        ax.bar(pos, values, width=bar_width, label=group)

    ax.set_xticks(x_pos + ((len(data) - 1) / 2) * bar_width)
    ax.set_xticklabels(xtick_labels)
    plt.xticks(rotation=90)

    ax.legend(legend)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)


def plot_dynamic(data, 
                 ylim=None, 
                 yscale="linear", 
                 xlabel="Time / hours", 
                 ylabel="Inoperability", 
                 figsize=(9, 5), 
                 dpi=300):
    """Helper function for plotting dynamic IIM data."""
    qt_data = data.head(-1)
    labels = qt_data.columns[1:]
    t_data = qt_data[qt_data.columns[0]].to_numpy()
    q_data = qt_data[qt_data.columns[1:]].to_numpy()

    _, ax = plt.subplots(figsize=figsize, dpi=dpi)
    colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(labels)))
    ax.set_prop_cycle("color", colors)

    for j in range(len(labels)):
        ax.plot(t_data, q_data[:, j], label=labels[j], linestyle="-", linewidth=2)

    if ylim:
        ax.set_ylim(ylim)
    plt.yscale(yscale)

    ax.yaxis.grid(color='gray', linestyle='dashed')

    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=5)

def plot_heatmap(df, vmin, vmax, xlabel=None, ylabel=None, cbar_label="Impact", dpi=300):
    """Helper function for creating heatmaps."""
    data = df.pivot("infra_i", "infra_j", "impact")
    _, ax = plt.subplots(dpi=dpi)
    sns.heatmap(data, linewidth=0.5, vmin=vmin, vmax=vmax, cbar_kws={"label": cbar_label}, ax=ax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


class PyDIIM:
    """Python wrapper for the DIIM code."""

    def __init__(self, user_config):
        self.config = {
            "job_name": "diim",
            "DIIM": {
                "amatrix_type": "input-output",
                "calc_mode": "demand",
                "amat_file": "diim_amat.csv",
                "kmat_file": "",
                "tau_file": "",
                "q0_file": "",
                "lambda_val": 0.01,
                "time_steps": 0 
            },
            "Perturbation": {
                "pinfra": [""],
                "cvalue": [0.0],
                "ptime": [[0, 0]] 
            }
        }
        self.__json_file = None
        if isinstance(user_config, str):
            self.__json_file = user_config
            self.config.update(self.__read_json_file())
        elif isinstance(user_config, dict):
            self.config.update(user_config)
            self.__json_file = self.config["job_name"] + ".json"
        self.__gen_json_file()

    def __read_json_file(self):
        """Read JSON input file."""
        data = {}
        with open(self.__json_file, "r") as f:
            data = json.load(f)
        return data

    def __gen_json_file(self):
        """Generate JSON input file."""
        with open(self.__json_file, "w") as f:
            json.dump(self.config, f)

    def run(self, run_type):
        """Run DIIM calculation."""
        try:
            output_file = self.config["job_name"] + ".csv"
            with open(output_file, "w") as f:
                subprocess.run(["diim_run", self.__json_file, run_type], stdout=f, check=True)
            return read_csv(output_file)
        except subprocess.CalledProcessError as exc:
            print(f"{exc}")
        except Exception as exc:
            print(f"{exc}")
