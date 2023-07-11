# Copyright (c) 2023 Stig Rune Sellevag
#
# This file is distributed under the MIT License. See the accompanying file
# LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
# and conditions.

"""Provides a Python wrapper for the DIIM C++ code."""

import numpy as np
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import AutoLocator, AutoMinorLocator


def score_to_interdependency(filename, score_scale):
    """Transform score values to interdependencies."""
    try:
        subprocess.run(["diim_gen", filename, str(score_scale)], check=True)
    except subprocess.FileNotFoundError as exc:
        print(f"{exc}")
    except subprocess.CalledProcessError as exc:
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

    def __init__(self, inp_file, config={}):
        self.__inp_file = inp_file
        self.config = {
            "amatrix_type": "input-output",
            "calc_mode": "demand",
            "amat_file": self.__inp_file.removesuffix(".inp") + ".csv",
            "kmat_file": None,
            "tau_file": None,
            "q0_file": None,
            "lambda_val": None,
            "time_steps": None,
            "pinfra": [],
            "cvalue": [],
            "ptime": [[]],
        }
        self.config.update(config)
        self.__gen_inp_file()

    def __gen_inp_file(self):
        """Generate input file."""
        with open(self.__inp_file, "w", encoding="latin1") as f:
            f.write("DIIM\n")
            f.write("  amatrix_type\n")
            f.write("    {0}\n".format(self.config["amatrix_type"]))
            f.write("  calc_mode\n")
            f.write("    {0}\n".format(self.config["calc_mode"]))
            f.write("  amat_file\n")
            f.write("    {0}\n".format(self.config["amat_file"]))
            if self.config["kmat_file"]:
                f.write("  kmat_file\n")
                f.write("    {0}\n".format(self.config["kmat_file"]))
            if self.config["tau_file"]:
                f.write("  tau_file\n")
                f.write("    {0}\n".format(self.config["tau_file"]))
            if self.config["q0_file"]:
                f.write("  q0_file\n")
                f.write("    {0}\n".format(self.config["q0_file"]))
            if self.config["lambda_val"]:
                f.write("  lambda\n")
                f.write("    {0}\n".format(self.config["lambda_val"]))
            if self.config["time_steps"]:
                f.write("  time_steps\n")
                f.write("    {0}\n".format(self.config["time_steps"]))
            f.write("End\n\n")

            if len(self.config["pinfra"]) > 0:
                f.write("Perturbation\n")
                f.write("  pinfra\n")
                f.write("    {0} [ ".format(len(self.config["pinfra"])))
                for i in self.config["pinfra"]:
                    f.write("{0} ".format(i))
                f.write("]\n")
                f.write("  cvalue\n")
                f.write("    {0} [ ".format(len(self.config["cvalue"])))
                for i in self.config["cvalue"]:
                    f.write("{0} ".format(i))
                f.write("]\n")
                f.write("  ptime\n")
                f.write("    {0}\n".format(len(self.config["ptime"])))
                for i in range(len(self.config["ptime"])):
                    f.write(
                        "    {0} {1}\n".format(
                            self.config["ptime"][i][0], self.config["ptime"][i][1]
                        )
                    )
                f.write("End\n")

    def run(self, run_type, output_file):
        """Run DIIM calculation."""
        try:
            with open(output_file, "w") as f:
                subprocess.run(["diim_run", self.__inp_file, run_type], stdout=f, check=True)
            return read_csv(output_file)
        except subprocess.FileNotFoundError as exc:
            print(f"{exc}")
        except subprocess.CalledProcessError as exc:
            print(f"{exc}")
        except Exception as exc:
            print(f"{exc}")
