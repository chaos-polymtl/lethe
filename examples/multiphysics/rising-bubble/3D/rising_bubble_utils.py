# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Functions for post-processing the 3D rising bubble problem (case #1)
"""
#############################################################################
#############################################################################
''' Importing Libraries '''

import os, re, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import MultipleLocator
import pyvista as pv
import importlib.util
import re

#############################################################################
''' Utility Functions '''
def setup_plot_params():
  '''
   Setup parameter values for plots
  '''
  plt.rcParams['figure.facecolor'] = 'white'
  plt.rcParams['figure.figsize'] = (10,8)
  plt.rcParams['lines.linewidth'] = 3
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
  plt.rcParams['xtick.top'] = True
  plt.rcParams['ytick.right'] = True
  plt.rcParams['xtick.direction'] = 'in'
  plt.rcParams['ytick.direction'] = 'in'
  plt.rcParams['font.size'] = '25'
  plt.rcParams['font.family']='DejaVu Serif'
  plt.rcParams['font.serif']='cm'
  plt.rcParams['savefig.bbox']='tight'

  plt.rcParams.update({
      'text.usetex': True,
      'text.latex.preamble': r'\usepackage{amsfonts}'
  })

def import_slice_extractor(lethe_path: str):
  '''
  Import the slice extractor to speed-up the extraction of bubble contours
  :param lethe_path: Local path to Lethe directory
  :return: Slice extractor module
  '''
  path_to_slice_extractor = f"{lethe_path}/contrib/postprocessing/extract-from-vtu/extract-slice-from-vtu.py"
  spec = importlib.util.spec_from_file_location("slice_extractor", path_to_slice_extractor)
  slice_extractor = importlib.util.module_from_spec(spec)
  sys.modules["slice_extractor"] = slice_extractor
  spec.loader.exec_module(slice_extractor)
  return slice_extractor

def base_format(path: str):
  '''
  Get base name of directory (the frequency value can be changed later)
  :param path: Path to directory
  :return: Directory name without the frequency
  '''
  if path is None:
    return None
  return re.sub(r'-f\d+', '-f{f:03d}', path)

def load_barycenter_series(folder: str):
  '''
  Load barycenter data
  :param folder:
  :return:
  '''
  fname = os.path.join(folder, "vof_barycenter_information.dat")
  t, x, y, z, vx, vy, vz = np.loadtxt(fname, skiprows=1, unpack=True)
  return {
    "t": t,
    "x": x,
    "y": y,
    "z": z,
    "vx": vx,
    "vy": vy,
    "vz": vz
  }

def load_method_data(has_method: bool, base_dir: str, f_values):
  '''
  Load data from ".dat" files
  :param has_method:
  :param base_dir:
  :param f_values:
  :return:
  '''
  bary, mass = {}, {}
  if not has_method:
    return bary, mass
  for i, f in enumerate(f_values):
    out = base_dir.format(f=f)
    mass[i] = os.path.join(out, "mass_conservation_information.dat")
    bary[i] = load_barycenter_series(out)
  return bary, mass

def frequency_labels(f_values):
  '''
  Frequency labels for legends
  :param f_values:
  :return:
  '''
  return [f"$f$=1/({f}" + r"$\Delta$t)" for f in f_values]

def extract_number(filename):
    match = re.search(r'\.(\d+)\.pvtu$', filename)
    return int(match.group(1)) if match else -1

def latest_sliced_contour(slice_extractor, out_dir: str, origin, normal, scalar="filtered_phase", iso=0.5):
  '''
  Pick latest .vtu (w/o "interface"/"sliced"), slice it, contour 0.5
  :param slice_extractor:
  :param out_dir:
  :param origin:
  :param normal:
  :param scalar:
  :param iso:
  :return:
  '''

  vtus = [f for f in os.listdir(out_dir) if (f.endswith('.pvtu') and "interface" not in f and "sliced" not in f)]
  latest = os.path.join(out_dir, max(vtus, key=extract_number))
  print(latest)

  slice_extractor.extract_slice(latest, origin, normal, out_dir)
  sliced = [os.path.join(out_dir, x) for x in os.listdir(out_dir) if "sliced" in x]
  last_sl = max(sliced, key=os.path.getctime)
  sim = pv.read(last_sl)
  sim.set_active_scalars(scalar)
  contours = sim.contour(np.array([iso]))
  pts = contours.points
  return pts[:, 0], pts[:, 1], pts[:, 2]


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


def read_refs():
  '''
  Read reference data from CSV files and save in a dictionary
  # [1] J. Adelsberger, P. Esser, M. Griebel, S. Groß, M. Klitz, and A. Rüttgers. 3D incompressible two-phase flow benchmark computations for rising droplets. 2014. Proceedings of the 11th World Congress on Computational Mechanics (WCCM XI), Barcelona, Spain, also available as INS Preprint No. 1401 and as IGPM Preprint No. 393.
  # [2] S. Groß, A. Reusken, Numerical methods for two-phase incompressible flows, Vol. 40 of Springer Series in Computational Mathematics, Springer, 2011
  # [3] Croce, R., Griebel, M. and Schweitzer, M.A. (2010), Numerical simulation of bubble and droplet deformation by a level set approach with surface tension in three dimensions. Int. J. Numer. Meth. Fluids, 62: 963-993. https://doi.org/10.1002/fld.2051
  # [4] OpenFOAM, The Open Source CFD Toolbox, User Guide Version 2.2.2, http://www.  openfoam.org (2013).
  :return: Dictionary of reference data values
  '''
  # FeatFlow [1]
  y_rise_velocity_featflow      = pd.read_csv("./reference-data/case_1/FeatFlow_results/FeatFlow_results/feat_flow_extracted_data/00_rise_velocity.csv")
  y_rise_velocity_featflow_zoom = pd.read_csv("./reference-data/case_1/FeatFlow_results/FeatFlow_results/feat_flow_extracted_data/00_rise_velocity_zoom.csv")
  sphericity_featflow           = pd.read_csv("./reference-data/case_1/FeatFlow_results/FeatFlow_results/feat_flow_extracted_data/02_bubble_sphericity.csv")
  sphericity_featflow_zoom      = pd.read_csv("./reference-data/case_1/FeatFlow_results/FeatFlow_results/feat_flow_extracted_data/02_bubble_sphericity_zoom.csv")
  # DROPS [2]
  y_rise_velocity_drops = pd.read_csv("./reference-data/case_1/ref_rise_velocity/TC1_riseVelocity_DROPS.dat", delim_whitespace=True, header=None);
  y_rise_velocity_drops.columns=["time","y_rise_velocity"]
  sphericity_drops      = pd.read_csv("./reference-data/case_1/ref_sphericity/TC1_sphericity_DROPS.dat", delim_whitespace=True, header=None);
  sphericity_drops.columns=["time","sphericity"]
  y_barycenter_drops    = pd.read_csv("./reference-data/case_1/ref_center_of_mass/TC1_centerOfMass_DROPS.dat", delim_whitespace=True, header=None);
  y_barycenter_drops.columns=["time","y_barycenter"]
  # NaSt3DGPF [3]
  y_rise_velocity_nast3dgpf = pd.read_csv("./reference-data/case_1/ref_rise_velocity/TC1_riseVelocity_NaSt3D.dat", delim_whitespace=True, header=None);
  y_rise_velocity_nast3dgpf.columns=["time","y_rise_velocity"]
  sphericity_nast3dgpf      = pd.read_csv("./reference-data/case_1/ref_sphericity/TC1_sphericity_NaSt3D.dat", delim_whitespace=True, header=None);
  sphericity_nast3dgpf.columns=["time","sphericity"]
  y_barycenter_nast3dgpf    = pd.read_csv("./reference-data/case_1/ref_center_of_mass/TC1_centerOfMass_NaSt3D.dat", delim_whitespace=True, header=None);
  y_barycenter_nast3dgpf.columns=["time","y_barycenter"]
  # OpenFOAM [4]
  y_rise_velocity_openfoam = pd.read_csv("./reference-data/case_1/ref_rise_velocity/TC1_riseVelocity_OpenFOAM.dat", delim_whitespace=True, header=None);
  y_rise_velocity_openfoam.columns=["time","y_rise_velocity"]
  sphericity_openfoam      = pd.read_csv("./reference-data/case_1/ref_sphericity/TC1_sphericity_OpenFOAM.dat", delim_whitespace=True, header=None);
  sphericity_openfoam.columns=["time","sphericity"]
  y_barycenter_openfoam    = pd.read_csv("./reference-data/case_1/ref_center_of_mass/TC1_centerOfMass_OpenFOAM.dat", delim_whitespace=True, header=None);
  y_barycenter_openfoam.columns=["time","y_barycenter"]
  return dict(
      featflow_velocity=y_rise_velocity_featflow, featflow_velocity_zoom=y_rise_velocity_featflow_zoom,
      featflow_sphericity=sphericity_featflow,     featflow_sphericity_zoom=sphericity_featflow_zoom,
      drops_velocity=y_rise_velocity_drops, drops_sphericity=sphericity_drops, drops_barycenter=y_barycenter_drops,
      nast3dgpf_velocity=y_rise_velocity_nast3dgpf, nast3dgpf_sphericity=sphericity_nast3dgpf, nast3dgpf_barycenter=y_barycenter_nast3dgpf,
      openfoam_velocity=y_rise_velocity_openfoam, openfoam_sphericity=sphericity_openfoam, openfoam_barycenter=y_barycenter_openfoam
  )

def save_tight(fig, path_base):
  '''

  :param fig:
  :param path_base:
  :return:
  '''
  fig.savefig(path_base + ".png", dpi=300, bbox_inches='tight')
  fig.savefig(path_base + ".svg", dpi=300, bbox_inches='tight')

#############################################################################
''' Legends '''
def make_legend_no_ref(save_path_base, f_values, c1_ref_linestyles, c1_alpha, lethe_linestyle, pde_colors,
                            geo_colors, pro_colors):
  '''

  :param save_path_base:
  :param f_values:
  :param c1_ref_linestyles:
  :param c1_alpha:
  :param lethe_linestyle:
  :param pde_colors:
  :param geo_colors:
  :param pro_colors:
  :return:
  '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  labels = frequency_labels(f_values)
  header = [plt.plot([], marker="", ls="")[0]]
  lethe_pde = [mlines.Line2D([], [], ls=lethe_linestyle, c=pde_colors[i]) for i in range(3)]
  lethe_geo = [mlines.Line2D([], [], ls=lethe_linestyle, c=geo_colors[i]) for i in range(3)]
  lethe_pro = [mlines.Line2D([], [], ls=lethe_linestyle, c=pro_colors[i]) for i in range(3)]

  lgd = ax.legend(
    handles=header + lethe_pde +
            header + lethe_geo +
            header + header + lethe_pro[1:],
    labels=["PDE-based"] + labels +
           ["Geometric"] + labels +
           ["Projection"] + [""] + labels[1:],
    ncol=3, loc="center", frameon=True, handler_map={tuple: HandlerTuple(ndivide=None)},
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

def make_legend_no_featflow(save_path_base, f_values, c1_ref_linestyles, c1_alpha, lethe_linestyle, pde_colors,
                            geo_colors, pro_colors):
  '''

  :param save_path_base:
  :param f_values:
  :param c1_ref_linestyles:
  :param c1_alpha:
  :param lethe_linestyle:
  :param pde_colors:
  :param geo_colors:
  :param pro_colors:
  :return:
  '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  labels = frequency_labels(f_values)
  header = [plt.plot([], marker="", ls="")[0]]
  lethe_pde = [mlines.Line2D([], [], ls=lethe_linestyle, c=pde_colors[i]) for i in range(3)]
  lethe_geo = [mlines.Line2D([], [], ls=lethe_linestyle, c=geo_colors[i]) for i in range(3)]
  lethe_pro = [mlines.Line2D([], [], ls=lethe_linestyle, c=pro_colors[i]) for i in range(3)]
  drops = mlines.Line2D([], [], ls=c1_ref_linestyles[0], c='k', alpha=c1_alpha)
  nast3d = mlines.Line2D([], [], ls=c1_ref_linestyles[2], c='k', alpha=c1_alpha)
  openfoam = mlines.Line2D([], [], ls=c1_ref_linestyles[3], c='k', alpha=c1_alpha)

  lgd = ax.legend(
    handles=header + lethe_pde +
            header + lethe_geo +
            header + header + lethe_pro[1:] +
            header + [drops, nast3d, openfoam],
    labels=["PDE-based"] + labels +
           ["Geometric"] + labels +
           ["Projection"] + [""] + labels[1:] +
           ["References"] + ["DROPS (2014)", "NaSt3DGPF (2014)", "OpenFOAM (2014)"],
    ncol=4, loc="center", frameon=True, handler_map={tuple: HandlerTuple(ndivide=None)},
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

def make_legend_with_featflow(save_path_base, f_values, c1_ref_linestyles, c1_ref_markers, c1_alpha, lethe_linestyle,
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
  labels = frequency_labels(f_values)
  header = [plt.plot([], marker="", ls="")[0]]

  lethe_pde = [mlines.Line2D([], [], ls=lethe_linestyle, c=pde_colors[i]) for i in range(3)]
  lethe_geo = [mlines.Line2D([], [], ls=lethe_linestyle, c=geo_colors[i]) for i in range(3)]
  lethe_pro = [mlines.Line2D([], [], ls=lethe_linestyle, c=pro_colors[i]) for i in range(3)]
  drops = mlines.Line2D([], [], ls=c1_ref_linestyles[0], c='k', alpha=c1_alpha, marker=c1_ref_markers[0])
  feat = mlines.Line2D([], [], ls=c1_ref_linestyles[1], c='k', alpha=c1_alpha, marker=c1_ref_markers[1])
  nast3d = mlines.Line2D([], [], ls=c1_ref_linestyles[2], c='k', alpha=c1_alpha, marker=c1_ref_markers[2])
  openfoam = mlines.Line2D([], [], ls=c1_ref_linestyles[3], c='k', alpha=c1_alpha, marker=c1_ref_markers[3])

  lgd = ax.legend(
    handles=header + lethe_pde + header +
            header + lethe_geo + header +
            header + header + lethe_pro[1:] + header +
            header + [drops, feat, nast3d, openfoam],
    labels=["PDE-based"] + labels + [""] +
           ["Geometric"] + labels + [""] +
           ["Projection"] + [""] + labels[1:] + [""] +
           ["References"] + ["DROPS (2014)", "FeatFlow (2019)", "NaSt3DGPF (2014)", "OpenFOAM (2014)"],
    ncol=4, loc="center", frameon=True, handler_map={tuple: HandlerTuple(ndivide=None)},
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


#############################################################################
''' Plots '''

def plot_barycenter(pde_bary, geo_bary, pro_bary, none_bary, refs, out_dir, f_values, case_number, line_width, c1_ref_linestyles, 
                    c1_color, c1_alpha, pde_colors, geo_colors, pro_colors, none_colors):
  '''
  
  :param pde_bary: 
  :param geo_bary: 
  :param pro_bary: 
  :param refs: 
  :param out_dir: 
  :param f_values: 
  :param case_number: 
  :param line_width: 
  :param c1_ref_linestyles: 
  :param c1_color: 
  :param c1_alpha: 
  :param pde_colors: 
  :param geo_colors: 
  :param pro_colors: 
  :return: 
  '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  fig_pro = plt.figure()
  ax_pro = fig_pro.add_subplot(111)

  axins = add_inset(ax, fig, (2.9, 3.005), (1.40, 1.48), locator=0.02, tick_right=True)
  axins_pro = add_inset(ax_pro, fig_pro, (2.9, 3.005), (1.35, 1.48), locator=0.05, tick_right=True)

  if pde_bary:
    for i, f in enumerate(f_values):
      ax.plot(pde_bary[i]["t"], pde_bary[i]["y"], '-', lw=line_width, color=pde_colors[i])
      axins.plot(pde_bary[i]["t"], pde_bary[i]["y"], '-', lw=line_width, color=pde_colors[i])
  if geo_bary:
    for i, f in enumerate(f_values):
      ax.plot(geo_bary[i]["t"], geo_bary[i]["y"], '-', lw=line_width, color=geo_colors[i])
      axins.plot(geo_bary[i]["t"], geo_bary[i]["y"], '-', lw=line_width, color=geo_colors[i])
  if pro_bary:
    for i, f in enumerate(f_values[1:]):
      ax_pro.plot(pro_bary[i + 1]["t"], pro_bary[i + 1]["y"], '-', lw=line_width, color=pro_colors[i + 1])
      axins_pro.plot(pro_bary[i + 1]["t"], pro_bary[i + 1]["y"], '-', lw=line_width, color=pro_colors[i + 1])
  if none_bary:
    i=0
    ax.plot(none_bary[i]["t"], none_bary[i]["y"], '-', lw=line_width, color=none_colors[i])
    axins.plot(none_bary[i]["t"], none_bary[i]["y"], '-', lw=line_width, color=none_colors[i])
    ax_pro.plot(none_bary[i]["t"], none_bary[i]["y"], '-', lw=line_width, color=none_colors[i])
    axins_pro.plot(none_bary[i]["t"], none_bary[i]["y"], '-', lw=line_width, color=none_colors[i])

  for a in (ax, ax_pro, axins, axins_pro):
    a.plot(refs["drops_barycenter"]["time"], refs["drops_barycenter"]["y_barycenter"], ls=c1_ref_linestyles[0], c=c1_color, lw=line_width,
           alpha=c1_alpha)
    a.plot(refs["nast3dgpf_barycenter"]["time"], refs["nast3dgpf_barycenter"]["y_barycenter"], ls=c1_ref_linestyles[2], c=c1_color, lw=line_width,
           alpha=c1_alpha)
    a.plot(refs["openfoam_barycenter"]["time"], refs["openfoam_barycenter"]["y_barycenter"], ls=c1_ref_linestyles[3], c=c1_color, lw=line_width,
           alpha=c1_alpha)

  ax.set_ylabel(r'Barycenter height [$\mathbb{L}$]')
  ax.set_xlabel(r'Rising time [$\mathbb{T}$]')
  ax_pro.set_ylabel(r'Barycenter height [$\mathbb{L}$]')
  ax_pro.set_xlabel(r'Rising time [$\mathbb{T}$]')
  save_tight(fig, f'{out_dir}/ymean-t-case{case_number}')
  save_tight(fig_pro, f'{out_dir}/ymean-t-case{case_number}-pro')
  plt.close(fig)
  plt.close(fig_pro)

def plot_rise_velocity(pde_bary, geo_bary, pro_bary, none_bary, refs, out_dir, f_values, case_number, line_width, c1_ref_linestyles,
                    c1_color, c1_alpha, pde_colors, geo_colors, pro_colors, none_colors):
  '''

  :param pde_bary:
  :param geo_bary:
  :param pro_bary:
  :param refs:
  :param out_dir:
  :param f_values:
  :param case_number:
  :param line_width:
  :param c1_ref_linestyles:
  :param c1_color:
  :param c1_alpha:
  :param pde_colors:
  :param geo_colors:
  :param pro_colors:
  :return:
  '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  fig_pro = plt.figure()
  ax_pro = fig_pro.add_subplot(111)
  axins = add_inset(ax, fig, (2.9, 3.0), (0.329, 0.355), box=(0.2, 0.08, 0.75, 0.75), loc="center", mark=(2, 4))

  if pde_bary:
    for i, f in enumerate(f_values):
      ax.plot(pde_bary[i]["t"], pde_bary[i]["vy"], '-', lw=line_width, color=pde_colors[i])
      axins.plot(pde_bary[i]["t"], pde_bary[i]["vy"], '-', lw=line_width, color=pde_colors[i])
  if geo_bary:
    for i, f in enumerate(f_values):
      ax.plot(geo_bary[i]["t"], geo_bary[i]["vy"], '-', lw=line_width, color=geo_colors[i])
      axins.plot(geo_bary[i]["t"], geo_bary[i]["vy"], '-', lw=line_width, color=geo_colors[i])
  if pro_bary:
    for i, f in enumerate(f_values[1:]):
      ax_pro.plot(pro_bary[i + 1]["t"], pro_bary[i + 1]["vy"], '-', lw=line_width, color=pro_colors[i + 1])
  if none_bary:
    i=0
    ax.plot(none_bary[i]["t"], none_bary[i]["vy"], '-', lw=line_width, color=none_colors[i])
    axins.plot(none_bary[i]["t"], none_bary[i]["vy"], '-', lw=line_width, color=none_colors[i])
    ax_pro.plot(none_bary[i]["t"], none_bary[i]["vy"], '-', lw=line_width, color=none_colors[i])


  for a in (ax, ax_pro, axins):
    a.plot(refs["drops_velocity"]["time"], refs["drops_velocity"]["y_rise_velocity"], ls=c1_ref_linestyles[0], c=c1_color,
           lw=line_width, alpha=c1_alpha)
    a.plot(refs["nast3dgpf_velocity"]["time"], refs["nast3dgpf_velocity"]["y_rise_velocity"], ls=c1_ref_linestyles[2], c=c1_color,
           lw=line_width, alpha=c1_alpha)
    a.plot(refs["openfoam_velocity"]["time"], refs["openfoam_velocity"]["y_rise_velocity"], ls=c1_ref_linestyles[3], c=c1_color,
           lw=line_width, alpha=c1_alpha)
  ax.plot(refs["featflow_velocity"]["t"], refs["featflow_velocity"]["rise_velocity"], ls=c1_ref_linestyles[1], c=c1_color,
          lw=line_width, alpha=c1_alpha)
  ax_pro.plot(refs["featflow_velocity"]["t"], refs["featflow_velocity"]["rise_velocity"], ls=c1_ref_linestyles[1], c=c1_color,
           lw=line_width, alpha=c1_alpha)
  axins.plot(refs["featflow_velocity"]["t"], refs["featflow_velocity"]["rise_velocity"], ls=c1_ref_linestyles[1], c=c1_color,
             lw=line_width, alpha=c1_alpha)

  ax.set_ylabel(r'Rise velocity [$\mathbb{LT}^{-1}$]');
  ax.set_xlabel(r'Rising time [$\mathbb{T}$]')
  ax_pro.set_ylabel(r'Rise velocity [$\mathbb{LT}^{-1}$]');
  ax_pro.set_xlabel(r'Rising time [$\mathbb{T}$]')
  save_tight(fig, f'{out_dir}/bubble-rise-velocity-case{case_number}')
  save_tight(fig_pro, f'{out_dir}/bubble-rise-velocity-case{case_number}-pro')
  plt.close(fig)
  plt.close(fig_pro)

def extract_contours(slice_extractor, has_pde, has_geo, has_pro, base_pde, base_geo, base_pro, f_values):
  '''
  Extract contours from
  Used in plot_contour_sets
  :param slice_extractor:
  :param has_pde:
  :param has_geo:
  :param has_pro:
  :param base_pde:
  :param base_geo:
  :param base_pro:
  :param f_values:
  :return:
  '''
  origin = np.array([0.5, 1, 0.5])
  normal = np.array([0, 0, 1])
  pde_c, geo_c, pro_c = {}, {}, {}
  if has_pde:
    for i, f in enumerate(f_values):
      x, y, z = latest_sliced_contour(slice_extractor, base_pde.format(f=f), origin, normal)
      pde_c[i] = {"x": x, "y": y, "z": z}
  if has_geo:
    for i, f in enumerate(f_values):
      x, y, z = latest_sliced_contour(slice_extractor, base_geo.format(f=f), origin, normal)
      geo_c[i] = {"x": x, "y": y, "z": z}
  if has_pro:
    for i, f in enumerate(f_values[1:]):
      x, y, z = latest_sliced_contour(slice_extractor, base_pro.format(f=f), origin, normal)
      pro_c[i] = {"x": x, "y": y, "z": z}
  return pde_c, geo_c, pro_c


def plot_contour_sets(pde_c, geo_c, pro_c, out_dir, f_values, case_number, line_width):
  '''

  :param pde_c:
  :param geo_c:
  :param pro_c:
  :param out_dir:
  :param f_values:
  :param case_number:
  :param line_width:
  :return:
  '''
  contours_info = []
  if pde_c:
    contours_info.append(("pde", pde_c, ['#54278f', '#756bb1', '#9e9ac8'], f_values))
  if geo_c:
    contours_info.append(("geo", geo_c, ['#a63603', '#e6550d', '#fd8d3c'], f_values))
  if pro_c:
    contours_info.append(("pro", pro_c, ['#31a354', '#74c476'], f_values[1:]))

  for name, data, cols, f_values in contours_info:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i, f in enumerate(f_values):
      d = data[i]
      ax.scatter(d["x"], d["y"], s=4, marker="s", color=cols[i], linewidths=line_width,
                 label=f'{"Projection" if name == "pro" else name.upper()} ($f$ = {f:03d})')
    ax.set_xlabel(r'x [$\mathbb{L}$]')
    ax.set_ylabel(r'y [$\mathbb{L}$]')
    ax.grid(which='major', color='grey', linestyle='--')
    ax.set_xlim([0.2, 0.8])
    ax.set_ylim([1.15, 1.75])
    save_tight(fig, f'{out_dir}/bubble-contour-{name}-case{case_number}')
    plt.close(fig)

def plot_global_mass_conservation(mass_pde, mass_geo, mass_pro, out_dir, f_values, case_number, line_width, pde_colors,
                                  geo_colors, pro_colors):
  '''

  :param mass_pde:
  :param mass_geo:
  :param mass_pro:
  :param out_dir:
  :param f_values:
  :param case_number:
  :param line_width:
  :param pde_colors:
  :param geo_colors:
  :param pro_colors:
  :return:
  '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  axins = add_inset(ax, fig, (2.9, 3.0), (0.94, 1.05), locator=0.05, tick_right=True)

  if mass_pde:
    for i, f in enumerate(f_values):
      df = pd.read_csv(mass_pde[i], header=0, sep=r'\s+')
      y = df['volume_fluid_1'] / df['volume_fluid_1'].iloc[0]
      ax.plot(df['time'], y, '-', linewidth=line_width, color=pde_colors[i])
      axins.plot(df['time'], y, '-', linewidth=line_width, color=pde_colors[i])
  if mass_geo:
    for i, f in enumerate(f_values):
      df = pd.read_csv(mass_geo[i], header=0, sep=r'\s+')
      y = df['volume_fluid_1'] / df['volume_fluid_1'].iloc[0]
      ax.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])
      axins.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])

  ax.set_ylabel(r'$V_\mathrm{phase\ fraction}/{V_\mathrm{0}}\ [-]$')
  ax.set_xlabel(r'Rising time [$\mathbb{T}$]')
  ax.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig, f'{out_dir}/global-mass-conservation-case{case_number}')
  plt.close(fig)

  if mass_pro:
    fig_pro = plt.figure()
    ax_pro = fig_pro.add_subplot(111)
    for i, f in enumerate(f_values[1:]):
      df = pd.read_csv(mass_pro[i + 1], header=0, sep=r'\s+')
      y = df['volume_fluid_1'] / df['volume_fluid_1'].iloc[0]
      ax_pro.plot(df['time'], y, '-', linewidth=line_width, color=pro_colors[i + 1])
    ax_pro.set_ylabel(r'$V_\mathrm{phase\ fraction}/{V_\mathrm{0}}\ [-]$')
    ax_pro.set_xlabel(r'Rising time [$\mathbb{T}$]')
    ax_pro.ticklabel_format(useOffset=False, style='plain', axis='y')
    save_tight(fig_pro, f'{out_dir}/global-mass-conservation-case{case_number}-pro')
    plt.close(fig_pro)


def plot_geometric_mass_conservation(mass_pde, mass_geo, mass_pro, mass_none, out_dir, f_values, case_number, line_width,
                                     pde_colors, geo_colors, pro_colors, none_colors):
  '''
  Plot geometric mass conservation
  :param mass_pde:
  :param mass_geo:
  :param mass_pro:
  :param out_dir:
  :param f_values:
  :param case_number:
  :param line_width:
  :param pde_colors:
  :param geo_colors:
  :param pro_colors:
  :return:
  '''
  fig = plt.figure()
  ax = fig.add_subplot(111)
  # axins = add_inset(ax, fig, (2.9,3.0), (0.965,1.02), locator=0.01, tick_right=True)

  if mass_pde:
      for i,f in enumerate(f_values):
          df = pd.read_csv(mass_pde[i], header=0, sep=r'\s+')
          y  = df['geometric_volume_fluid_1']/df['geometric_volume_fluid_1'].iloc[0]
          ax.plot(df['time'], y, '-', linewidth=line_width, color=pde_colors[i])
          # axins.plot(df['time'], y, '-', linewidth=line_width, color=pde_colors[i])
  if mass_geo:
      for i,f in enumerate(f_values):
          df = pd.read_csv(mass_geo[i], header=0, sep=r'\s+')
          y  = df['geometric_volume_fluid_1']/df['geometric_volume_fluid_1'].iloc[0]
          ax.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])
          # axins.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])
  if mass_none:
    i=0
    df = pd.read_csv(mass_none[i], header=0, sep=r'\s+')
    y  = df['geometric_volume_fluid_1']/df['geometric_volume_fluid_1'].iloc[0]
    ax.plot(df['time'], y, '-', linewidth=line_width, color=none_colors[i])
    # axins.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])

  ax.set_ylabel(r'$V/{V_\mathrm{0}}\ [-]$')
  ax.set_xlabel(r'Rising time [$\mathbb{T}$]')
  ax.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig, f'{out_dir}/geo-mass-conservation-case{case_number}')
  plt.close(fig)

  if mass_pro:
      fig_pro = plt.figure()
      ax_pro = fig_pro.add_subplot(111)
      for i,f in enumerate(f_values[1:]):
          df = pd.read_csv(mass_pro[i+1], header=0, sep=r'\s+')
          y  = df['geometric_volume_fluid_1']/df['geometric_volume_fluid_1'].iloc[0]
          ax_pro.plot(df['time'], y, '-', linewidth=line_width, color=pro_colors[i+1])
  if mass_none:
    i=0
    df = pd.read_csv(mass_none[i], header=0, sep=r'\s+')
    y  = df['geometric_volume_fluid_1']/df['geometric_volume_fluid_1'].iloc[0]
    ax_pro.plot(df['time'], y, '-', linewidth=line_width, color=none_colors[i])
    # axins.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])

  ax_pro.set_ylabel(r'$V/{V_\mathrm{0}}\ [-]$')
  ax_pro.set_xlabel(r'Rising time [$\mathbb{T}$]')
  ax_pro.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig_pro, f'{out_dir}/geo-mass-conservation-case{case_number}-pro')
  plt.close(fig_pro)

def plot_sphericity(mass_pde, mass_geo, mass_pro, mass_none, refs, out_dir, f_values, case_number, line_width, c1_ref_linestyles,
                    c1_color, c1_alpha, pde_colors, geo_colors, pro_colors, none_colors):
  '''

  :param mass_pde:
  :param mass_geo:
  :param mass_pro:
  :param refs:
  :param out_dir:
  :param f_values:
  :param case_number:
  :param line_width:
  :param c1_ref_linestyles:
  :param c1_color:
  :param c1_alpha:
  :param pde_colors:
  :param geo_colors:
  :param pro_colors:
  :return:
  '''

  fig = plt.figure()
  ax = fig.add_subplot(111)
  fig_pro = plt.figure()
  ax_pro = fig_pro.add_subplot(111)
  axins = add_inset(ax, fig, (2.5, 3.0), (0.957, 0.965),
                    box=(0.4, 0.32, 0.55, 0.55), loc="center", mark=(3, 1), locator=0.005)
  axins.xaxis.tick_top()

  if mass_pde:
    for i, f in enumerate(f_values):
      df = pd.read_csv(mass_pde[i], header=0, sep=r'\s+')
      area_eq = np.pi ** (1 / 3) * (6 * df["geometric_volume_fluid_1"]) ** (2 / 3)
      y = area_eq / df["surface_fluid_1"]
      ax.plot(df['time'], y, '-', linewidth=line_width, color=pde_colors[i])
      axins.plot(df['time'], y, '-', linewidth=line_width, color=pde_colors[i])
  if mass_geo:
    for i, f in enumerate(f_values):
      df = pd.read_csv(mass_geo[i], header=0, sep=r'\s+')
      area_eq = np.pi ** (1 / 3) * (6 * df["geometric_volume_fluid_1"]) ** (2 / 3)
      y = area_eq / df["surface_fluid_1"]
      ax.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])
      axins.plot(df['time'], y, '-', linewidth=line_width, color=geo_colors[i])
  if mass_pro:
    for i, f in enumerate(f_values[1:]):
      df = pd.read_csv(mass_pro[i + 1], header=0, sep=r'\s+')
      area_eq = np.pi ** (1 / 3) * (6 * df["geometric_volume_fluid_1"]) ** (2 / 3)
      y = area_eq / df["surface_fluid_1"]
      ax_pro.plot(df['time'], y, '-', linewidth=line_width, color=pro_colors[i + 1])
  if mass_none:
    i=0
    df = pd.read_csv(mass_none[i], header=0, sep=r'\s+')
    area_eq = np.pi ** (1 / 3) * (6 * df["geometric_volume_fluid_1"]) ** (2 / 3)
    y = area_eq / df["surface_fluid_1"]
    ax.plot(df['time'], y, '-', linewidth=line_width, color=none_colors[i])
    axins.plot(df['time'], y, '-', linewidth=line_width, color=none_colors[i])
    ax_pro.plot(df['time'], y, '-', linewidth=line_width, color=none_colors[i])


  for a in (ax, ax_pro, axins):
    a.plot(refs["drops_sphericity"]["time"], refs["drops_sphericity"]["sphericity"], ls=c1_ref_linestyles[0], c=c1_color,
           lw=line_width, alpha=c1_alpha)
    a.plot(refs["nast3dgpf_sphericity"]["time"], refs["nast3dgpf_sphericity"]["sphericity"], ls=c1_ref_linestyles[2], c=c1_color,
           lw=line_width, alpha=c1_alpha)
    a.plot(refs["openfoam_sphericity"]["time"], refs["openfoam_sphericity"]["sphericity"], ls=c1_ref_linestyles[3], c=c1_color, lw=line_width,
           alpha=c1_alpha)
    a.plot(refs["featflow_sphericity"]["t"], refs["featflow_sphericity"]["sphericity"], ls=c1_ref_linestyles[1], c=c1_color,
            lw=line_width, alpha=c1_alpha)

  ax.set_ylabel(r'$A_\mathrm{eq}/{A_\Gamma}\ [-]$')
  ax.set_xlabel(r'Rising time [$\mathbb{T}$]')
  ax.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig, f'{out_dir}/bubble-sphericity-case{case_number}')
  plt.close(fig)

  ax_pro.set_ylabel(r'$A_\mathrm{eq}/{A_\Gamma}\ [-]$')
  ax_pro.set_xlabel(r'Rising time [$\mathbb{T}$]')
  ax_pro.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig_pro, f'{out_dir}/bubble-sphericity-case{case_number}-pro')
  plt.close(fig_pro)


def load_cfl_evolution(has_method: bool, base_dir: str, f_values):
  '''

  :param has_method:
  :param base_dir:
  :param f_values:
  :return:
  '''
  cfl_evo = {}
  if not has_method:
      return
  for i, f in enumerate(f_values):
    out = base_dir.format(f=f)
    cfl_evo[i] = os.path.join(out, "cfl_evolution.csv")
  return cfl_evo

def plot_cfl_evolution(cfl_pde, cfl_geo, cfl_pro, out_dir, f_values, case_number, line_width, pde_colors, geo_colors, pro_colors):
  fig = plt.figure()
  ax = fig.add_subplot(111)
  fig_pro = plt.figure()
  ax_pro = fig_pro.add_subplot(111)
  
  if (cfl_pde):
    for i, f in enumerate(f_values):
      df = pd.read_csv(cfl_pde[i])
      ax.plot(df['Time'], df['CFL'],'-', linewidth=line_width, color=pde_colors[i])
  if (cfl_geo):
    for i, f in enumerate(f_values):
      df = pd.read_csv(cfl_geo[i])
      ax.plot(df['Time'], df['CFL'],'-', linewidth=line_width, color=geo_colors[i])
  if (cfl_pro):
    for i, f in enumerate(f_values):
      df = pd.read_csv(cfl_pro[i])
      ax_pro.plot(df['Time'], df['CFL'],'-', linewidth=line_width, color=pro_colors[i])

  ax.set_ylabel(r"CFL $\left(\frac{u\Delta t}{h}\right)$ [-]")
  ax.set_xlabel(r"Time [$\mathbb{T}$]")
  ax.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig, f'{out_dir}/cfl-evolution-case{case_number}')
  plt.close(fig)

  ax_pro.set_ylabel(r"CFL $\left(\frac{u\Delta t}{h}\right)$ [-]")
  ax_pro.set_xlabel(r"Time [$\mathbb{T}$]")
  ax_pro.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig_pro, f'{out_dir}/cfl-evolution-case{case_number}-pro')
  plt.close(fig_pro)
