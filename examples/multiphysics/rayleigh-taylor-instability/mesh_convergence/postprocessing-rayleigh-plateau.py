# SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for Rayleigh-Plateau example

"""
#############################################################################

'''Importing Libraries'''
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import MultipleLocator

import pyvista as pv
import argparse
import os
import pandas as pd
#############################################################################

'''Plot formating'''

from cycler import cycler

colors=[]
pro_colors = ['#006d2c','#31a354','#74c476','#bae4b3']
pde_colors = ['#54278f','#756bb1','#9e9ac8','#cbc9e2']
geo_colors = ['#a63603','#e6550d','#fd8d3c','#fdbe85']
none_colors = ['#d7b5d8','#df65b0','#dd1c77','#980043']


plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.markersize'] = '11'
plt.rcParams['markers.fillstyle'] = "none"
plt.rcParams['lines.markeredgewidth'] = 2
plt.rcParams['legend.columnspacing'] = 2
plt.rcParams['legend.handlelength'] = 2.8
plt.rcParams['legend.handletextpad'] = 0.2
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.fontsize'] = '18'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['font.size'] = '25'
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['savefig.bbox']='tight'

plt.rcParams.update({
    'text.usetex': True,
    'text.latex.preamble': r'\usepackage{amsfonts}'
})

#############################################################################
def save_tight(fig, path_base):
  '''

  :param fig:
  :param path_base:
  :return:
  '''
  fig.savefig(path_base + ".png", dpi=300, bbox_inches='tight')
  fig.savefig(path_base + ".svg", dpi=300, bbox_inches='tight')

def make_legend(save_path_base, lethe_linestyle,
                              pde_colors, geo_colors, pro_colors):
  '''

  :param save_path_base:
  :param f_values:
  :param c1_ref_linestyles:
  :param c1_ref_markers:
  :param c1_alpha:
  :param lethe_linestyle:
  :param pde_colors:
  :param geo_colors:
  :param pro_colors:
  :return:
  '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  labels = ["Coarse", "Medium", "Fine", "Extra-fine"]
  header = [plt.plot([], marker="", ls="")[0]]

  lethe_pde = [mlines.Line2D([], [], ls=lethe_linestyle, c=pde_colors[i]) for i in range(4)]
  lethe_geo = [mlines.Line2D([], [], ls=lethe_linestyle, c=geo_colors[i]) for i in range(4)]

  lgd = ax.legend(
    handles=header + lethe_pde +
            header + lethe_geo,
    labels=["PDE-based"] + labels+
           ["Geometric"] + labels,
    ncol=2, loc="center", frameon=True, handler_map={tuple: HandlerTuple(ndivide=None)},
    handletextpad=0.2)
  for vpack in lgd._legend_handle_box.get_children():
    for hpack in vpack.get_children()[:1]:
      hpack.get_children()[0].set_width(0)
  ax.axis('off')
  fig.canvas.draw()
  bbox = lgd.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
  fig.savefig(save_path_base + ".png", dpi=300, bbox_inches=bbox)
  fig.savefig(save_path_base + ".svg", dpi=300, bbox_inches=bbox)
  plt.close(fig)

def add_inset(ax, fig, xlims, ylims, locator=None, box=(0.8, 0.1, 0.8, 0.8), loc="center",
              mark=(2,3), tick_right=False, tick_top=False):
  '''

  :param ax:
  :param fig:
  :param xlims:
  :param ylims:
  :param locator:
  :param box:
  :param loc:
  :param mark:
  :param tick_right:
  :param tick_top:
  :return:
  '''
  axins = inset_axes(ax, width="60%", height="60%",
                     bbox_to_anchor=box, bbox_transform=fig.transFigure, loc=loc)
  axins.set_xlim(*xlims); axins.set_ylim(*ylims)
  if locator:
    axins.yaxis.set_major_locator(MultipleLocator(locator))
  mark_inset(ax, axins, loc1=mark[0], loc2=mark[1], fc="none", ec="0.45")
  if tick_right:
    axins.yaxis.tick_right()
  if tick_top:
    axins.xaxis.tick_top()
  return axins

parser = argparse.ArgumentParser(description='Arguments for the postprocess of the Rayleigh-Plateau case')
parser.add_argument("-f", "--file", type=str, nargs='+', help="Path to the output folder for the results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=True)
args, leftovers=parser.parse_known_args()

make_legend(f'./legend',  '-', pde_colors, geo_colors, pro_colors)

for k in range(len(args.file)):
    fig_volume = plt.figure()
    ax_volume = fig_volume.add_subplot(111)

    fig_surface = plt.figure()
    ax_surface = fig_surface.add_subplot(111)

    output_dir_base = args.file[k]

    ext = ""
    if "geo" in output_dir_base:
        colors = geo_colors
        ext = "geo"
    elif "pro" in output_dir_base:
        colors = pro_colors
        ext = "pro"
    elif "pde" in output_dir_base:
        colors = pde_colors
        ext = "pde"
    else:
        colors = none_colors
        ext = "none"


    plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
    j=0
    for i in range(6,10):
        output_dir = output_dir_base.replace("ref_max05", "ref_max0"+str(i))
        filename_mass = output_dir + "/mass_conservation_information.dat"

        df_mass = pd.read_csv(filename_mass, header=0, sep=r'\s+')

        label = "Ref. max " + str(i)
        ax_volume.plot(df_mass['time'], 
            df_mass['geometric_volume_fluid_1'] / df_mass['geometric_volume_fluid_1'].iloc[0],
            "-", label=label, linewidth=2,color=colors[j])

        ax_surface.plot(df_mass['time'], df_mass['surface_fluid_1'] / df_mass['surface_fluid_1'].iloc[0],
            "-", label=label, linewidth=2,color=colors[j])
        j += 1

    ax_volume.set_ylabel(r'$V/{V_\mathrm{0}}[-]$')
    ax_volume.set_xlabel(r'Time [s]')
    ax_volume.set_xlim([-0.002, 0.08])
    ax_volume.set_ylim([0.78, 1.08])

    ax_surface.set_ylabel(r'$A_\Gamma/A_{\Gamma,0}[-]$')
    ax_surface.set_xlabel(r'Time [s]')
    ax_surface.set_xlim([-0.002, 0.08])
    ax_surface.set_ylim([0.76, 1.02])

    fig_volume_filename = './' + ext + '-volume-evolution.png'
    fig_volume.savefig(fig_volume_filename,dpi=300, bbox_inches='tight')

    fig_surface_filename = './' + ext + '-surface-evolution.png'

    fig_surface.savefig(fig_surface_filename,dpi=300, bbox_inches='tight')

if (len(args.file) > 1):
  ext = "comparison"
  fig_volume = plt.figure()
  ax_volume = fig_volume.add_subplot(111)

  fig_surface = plt.figure()
  ax_surface = fig_surface.add_subplot(111)
  for k in range(len(args.file)):

    output_dir_base = args.file[k]

    i=9
    if "geo" in output_dir_base:
        colors = geo_colors
    elif "pro" in output_dir_base:
        colors = pro_colors
    elif "pde" in output_dir_base:
        colors = pde_colors
    else:
        colors = none_colors

    plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
    j=3

    output_dir = output_dir_base.replace("ref_max05", "ref_max0"+str(i))
    filename_mass = output_dir + "/mass_conservation_information.dat"

    df_mass = pd.read_csv(filename_mass, header=0, sep=r'\s+')

    label = "Ref. max " + str(i)
    ax_volume.plot(df_mass['time'], 
        df_mass['geometric_volume_fluid_1'] / df_mass['geometric_volume_fluid_1'].iloc[0],
        "-", label=label, linewidth=2,color=colors[j])

    ax_surface.plot(df_mass['time'], df_mass['surface_fluid_1'] / df_mass['surface_fluid_1'].iloc[0],
        "-", label=label, linewidth=2,color=colors[j])
    j += 1

  ax_volume.set_ylabel(r'$V/{V_\mathrm{0}}[-]$')
  ax_volume.set_xlabel(r'Time [s]')
  ax_volume.set_xlim([-0.002, 0.08])
  ax_volume.set_ylim([0.78, 1.08])

  ax_surface.set_ylabel(r'$A_\Gamma/A_{\Gamma,0}[-]$')
  ax_surface.set_xlabel(r'Time [s]')
  ax_surface.set_xlim([-0.002, 0.08])
  ax_surface.set_ylim([0.76, 1.02])

  fig_volume_filename = './' + ext + '-volume-evolution.png'
  fig_volume.savefig(fig_volume_filename,dpi=300, bbox_inches='tight')

  fig_surface_filename = './' + ext + '-surface-evolution.png'

  fig_surface.savefig(fig_surface_filename,dpi=300, bbox_inches='tight')


