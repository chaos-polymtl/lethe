# SPDX-FileCopyrightText: Copyright (c) 2022-2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for rising-bubble example
"""
#############################################################################

'''Importing Libraries'''
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import pyvista as pv
import argparse
import os
import pandas as pd
import importlib.util
import re
#############################################################################
def save_tight(fig, path_base):
  '''

  :param fig:
  :param path_base:
  :return:
  '''
  fig.savefig(path_base + ".png", dpi=300, bbox_inches='tight')
  fig.savefig(path_base + ".svg", dpi=300, bbox_inches='tight')

def base_format(path: str):
  '''
  Get base name of directory (the frequency value can be changed later)
  :param path: Path to directory
  :return: Directory name without the frequency
  '''
  if path is None:
    return None
  return re.sub(r'-f\d+', '-f{f:03d}', path)
  
def extract_number(filename):
    match = re.search(r'\.(\d+)\.pvtu$', filename)
    return int(match.group(1)) if match else -1

def frequency_labels(f_values):
  '''
  Frequency labels for legends
  :param f_values:
  :return:
  '''
  return [f"$f$=1/({f}" + r"$\Delta$t)" for f in f_values]

def make_legend(save_path_base, f_values, c1_ref_linestyles, c1_alpha, lethe_linestyle, pde_colors,
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
  analytic = mlines.Line2D([], [], ls=c1_ref_linestyles[0], c='k', alpha=c1_alpha)
  empty_lgd = mlines.Line2D([], [], ls=c1_ref_linestyles[0], c='k', alpha=0)

  lgd = ax.legend(
    handles=header + lethe_pde +
            header + lethe_geo +
            header + header + lethe_pro[1:] +
            header + [analytic,empty_lgd,empty_lgd],
    labels=["PDE-based"] + labels +
           ["Geometric"] + labels +
           ["Projection"] + [""] + labels[1:] +
           ["Reference"] + ["Analytical solution", " ", " "],
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

def plot_cfl_evolution(cfl_pde, cfl_geo, cfl_pro, out_dir, f_values, line_width, pde_colors, geo_colors, pro_colors):
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
  ax.set_xlabel(r"Time [T]")
  ax.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig, f'{out_dir}/cfl-evolution')
  plt.close(fig)

  ax_pro.set_ylabel(r"CFL $\left(\frac{u\Delta t}{h}\right)$ [-]")
  ax_pro.set_xlabel(r"Time [T]")
  ax_pro.ticklabel_format(useOffset=False, style='plain', axis='y')
  save_tight(fig_pro, f'{out_dir}/cfl-evolution-pro')
  plt.close(fig_pro)

'''Plot formating'''

from cycler import cycler

colors=['#008c66','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02']
pro_colors = ['#006d2c','#31a354','#74c476','#bae4b3']
pde_colors = ['#54278f','#756bb1','#9e9ac8','#cbc9e2']
geo_colors = ['#a63603','#e6550d','#fd8d3c','#fdbe85']

c1_ref_linestyles = ['--']
c1_ref_markers = [0]
c1_ref_label = ["Analytic solution"]
c1_alpha = 0.8
c1_color = 'k'
lethe_label = "Lethe (2025)"
lethe_linestyle = "-"

plt.rcParams['axes.prop_cycle'] = cycler(color = colors)
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.figsize'] = (10,8)
plt.rcParams['lines.linewidth'] = 3
plt.rcParams['lines.markersize'] = '11'
# plt.rcParams['markers.fillstyle'] = "none"
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

#############################################################################


parser = argparse.ArgumentParser(description='Arguments for the validation of the 3D capillary migration case')
parser.add_argument("-p", "--pro", type=str, help="Path to the output folder for the projection results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=False)
parser.add_argument("-g", "--geo", type=str, help="Path to the output folder for the geometric results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=False)
parser.add_argument("-a", "--pde", type=str, help="Path to the output folder for the algebraic results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=False)
parser.add_argument("-l", "--lethe_path", type=str, help="Path to the Lethe's main folder", required=True)
parser.add_argument("-sf", "--save_figure_path", type=str, help="Path to where figure should be saved", required=False, default="./")
args, leftovers=parser.parse_known_args()
parser.add_argument("-fs", "--frequency_sweep", type=bool, help="Post-process frequency sweep if set to true", required=False)
args, leftovers=parser.parse_known_args()

# Create save folder if it does not exist
os.makedirs(args.save_figure_path, exist_ok=True)

# Import VTU slicer
path_to_slice_extractor=f"{args.lethe_path}/contrib/postprocessing/extract-from-vtu/extract-slice-from-vtu.py"
spec = importlib.util.spec_from_file_location("slice_extractor", path_to_slice_extractor)
slice_extractor = importlib.util.module_from_spec(spec)
sys.modules["slice_extractor"] = slice_extractor
spec.loader.exec_module(slice_extractor)

has_pro = args.pro is not None
has_geo  = args.geo  is not None
has_pde = args.pde is not None

if (args.frequency_sweep):
  f_values = [ 1, 10, 100]
  base_output_dir_pde = re.sub(r'-f\d+', '-f{f:03d}', args.pde)
  base_output_dir_geo = re.sub(r'-f\d+', '-f{f:03d}', args.geo)
  base_output_dir_pro = re.sub(r'-f\d+', '-f{f:03d}', args.pro)

  # Store data by frequency value
  pde_barycenter_data = {}
  pde_sphericity_data = {}
  pde_mass_geo_data = {}
  pde_contour_data = {}
  geo_barycenter_data = {}
  geo_sphericity_data = {}
  geo_mass_geo_data = {}
  geo_contour_data = {}
  pro_barycenter_data = {}
  pro_sphericity_data = {}
  pro_mass_geo_data = {}
  pro_contour_data = {}


make_legend(f'{args.save_figure_path}/legend', f_values,  c1_ref_linestyles,
                          c1_alpha, lethe_linestyle, pde_colors, geo_colors, pro_colors)

# Number of datasets
n = 0

if has_pro:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      output_dir_pro = base_output_dir_pro.format(f=f)
      filename_barycenter_pro = output_dir_pro + "/barycenter_information.dat"
      pro_mass_geo_data[i] = output_dir_pro + "/mass_conservation_information.dat"
      t_pro, x_pro, y_pro, z_pro, vx_pro, vy_pro, vz_pro = np.loadtxt(filename_barycenter_pro, skiprows=1, unpack=True)

      # Store in dictionary
      pro_barycenter_data[i] = {
        "t": t_pro,
        "x": x_pro,
        "y": y_pro,
        "z": z_pro,
        "vx": vx_pro,
        "vy": vy_pro,
        "vz": vz_pro
      }
      n += 1
  else:
    output_dir_pro = args.pro
    filename_barycenter_pro = output_dir_pro + "/barycenter_information.dat"
    filename_mass_pro = output_dir_pro + "/mass_conservation_information.dat"
    t_pro, x_pro, y_pro,z_pro, vx_pro, vy_pro , vz_pro = np.loadtxt(filename_barycenter_pro, skiprows=1, unpack=True)
    n += 1
if has_geo:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      output_dir_geo = base_output_dir_geo.format(f=f)
      filename_barycenter_geo = output_dir_geo + "/barycenter_information.dat"
      geo_mass_geo_data[i] = output_dir_geo + "/mass_conservation_information.dat"
      t_geo, x_geo, y_geo, z_geo, vx_geo, vy_geo, vz_geo = np.loadtxt(filename_barycenter_geo, skiprows=1, unpack=True)

      # Store in dictionary
      geo_barycenter_data[i] = {
        "t": t_geo,
        "x": x_geo,
        "y": y_geo,
        "z": z_geo,
        "vx": vx_geo,
        "vy": vy_geo,
        "vz": vz_geo
      }
      n += 1
  else:
    output_dir_geo = args.geo
    filename_barycenter_geo = output_dir_geo + "/barycenter_information.dat"
    filename_mass_geo = output_dir_geo + "/mass_conservation_information.dat"
    t_geo, x_geo, y_geo, z_geo, vx_geo, vy_geo, vz_geo = np.loadtxt(filename_barycenter_geo, skiprows=1, unpack=True)
    n += 1

if has_pde:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      output_dir_pde = base_output_dir_pde.format(f=f)
      filename_barycenter_pde = output_dir_pde + "/barycenter_information.dat"
      pde_mass_geo_data[i] = output_dir_pde + "/mass_conservation_information.dat"
      t_pde, x_pde, y_pde, z_pde, vx_pde, vy_pde, vz_pde = np.loadtxt(filename_barycenter_pde, skiprows=1, unpack=True)

      # Store in dictionary
      pde_barycenter_data[i] = {
        "t": t_pde,
        "x": x_pde,
        "y": y_pde,
        "z": z_pde,
        "vx": vx_pde,
        "vy": vy_pde,
        "vz": vz_pde
      }
      n += 1
  else:
    output_dir_pde = args.pde
    filename_barycenter_pde = output_dir_pde + "/barycenter_information.dat"
    filename_mass_pde = output_dir_pde + "/mass_conservation_information.dat"
    t_pde, x_pde, y_pde,  z_pde, vx_pde, vy_pde, vz_pde = np.loadtxt(filename_barycenter_pde, skiprows=1, unpack=True)
    n += 1

"""
 Analytic solution 
"""

radius = 0.25

# dgamma/dT*|grad T|
gamma_prime = -1

lambda_ratio = 1
beta_ratio = 1

rho = 1

nu_ext = 1
nu_int = nu_ext*lambda_ratio

mu_ext = nu_ext*rho
mu_int = nu_int*rho

t_analytical = np.linspace(0.0,3.0, 25)

vx_analytical = -2*radius*gamma_prime/(mu_ext*(2+3*lambda_ratio)*(2+beta_ratio))*np.ones(len(t_analytical))

print(vx_analytical)
# """
#  Construct legends 
# """

# ############################
# # No Featflow
# ############################
# fig00 = plt.figure()
# ax00 = fig00.add_subplot(111)

# legend = []
# labels = []
# n_columns = 3

# header = [plt.plot([],marker="", ls="")[0]]

# for i, f in enumerate(f_values):
#   row = [] # clear row

#   pde = mlines.Line2D([], [], color=pde_colors[i])
#   geo = mlines.Line2D([], [], color=geo_colors[i])
#   pro = mlines.Line2D([], [], color=pro_colors[i])
#   row = (pde, geo, pro)

#   legend.append(row)
#   labels.append(f"$f$={1/f:<1.2f}")

# # Lethe
# lethe_pde_lgd = mlines.Line2D([], [], ls=lethe_linestyle, c=pde_colors[1])
# lethe_geo_lgd = mlines.Line2D([], [], ls=lethe_linestyle, c=geo_colors[1])
# lethe_pro_lgd = mlines.Line2D([], [], ls=lethe_linestyle, c=pro_colors[1])

# # Reference solutions
# analytic_lgd = mlines.Line2D([], [], ls=c1_ref_linestyles[0], c='k', alpha=c1_alpha)
# empty_lgd = mlines.Line2D([], [], ls=c1_ref_linestyles[0], c='k', alpha=0)

# leg_no_featflow = ax00.legend(handles= header + legend +
#                                       header + [lethe_pde_lgd,lethe_geo_lgd,lethe_pro_lgd]+
#                                       header + [analytic_lgd,empty_lgd,empty_lgd],
#                              labels= [f"Reinitialization \n frequencies $(f)$"] + labels +
#                                      [lethe_label] + ["PDE-based", "Geometric", "Projection"]+
#                                      ["Reference"] + [c1_ref_label[0]," "," "],
#                              ncol=n_columns,
#                              loc="center",
#                              frameon=True,
#                              handler_map={tuple: HandlerTuple(ndivide=None)},
#                              handletextpad=0.2)

# for vpack in leg_no_featflow._legend_handle_box.get_children():
#   for hpack in vpack.get_children()[:1]:
#     hpack.get_children()[0].set_width(0)

# ax00.axis('off')

# fig00.canvas.draw()
# bbox = leg_no_featflow.get_window_extent()
# bbox = bbox.transformed(fig00.dpi_scale_trans.inverted())

# fig00.savefig(f'{args.save_figure_path}/legend.png',
#               dpi=300, bbox_inches=bbox)
# fig00.savefig(f'{args.save_figure_path}/legend.svg',
#               dpi=300, bbox_inches=bbox)   
                         
#************************************************************
# Plot the migation velocity
#************************************************************

fig0 = plt.figure()
ax0 = fig0.add_subplot(111)

fig0_proj = plt.figure()
ax0_proj = fig0_proj.add_subplot(111)

line_width = 3.0
if n == 1:
    line_width = 3.0

if has_pde:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      ax0.plot(pde_barycenter_data[i]["t"], pde_barycenter_data[i]["vx"], '-', lw=line_width, label=f"PDE-based ($f$={f:03d})", color = pde_colors[i])
  else:
    ax0.plot(t_pde, y_pde, '-', lw=line_width, label="PDE-based")
if has_geo:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      ax0.plot(geo_barycenter_data[i]["t"], geo_barycenter_data[i]["vx"], '-', lw=line_width, label=f"Geometric ($f$={f:03d})", color = geo_colors[i])
  else:
    ax0.plot(t_geo, y_geo, '-', lw=line_width, label="Geometric")
if has_pro:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      if (f == 1):
        continue
      ax0_proj.plot(pro_barycenter_data[i]["t"], pro_barycenter_data[i]["vx"], '-', lw=line_width, label=f"Projection ($f$={f:03d})", color = pro_colors[i])
  else:
    ax0_proj.plot(t_pro, y_pro, '-', lw=line_width, label="Projection")

ax0.plot(t_analytical, vx_analytical, ls=c1_ref_linestyles[0], c=c1_color, lw=line_width, alpha=c1_alpha, label="Analytic")
ax0_proj.plot(t_analytical, vx_analytical, ls=c1_ref_linestyles[0], c=c1_color, lw=line_width, alpha=c1_alpha, label="Analytic")

# inset axes
axins = inset_axes(
  ax0, width="60%", height="60%",
  bbox_to_anchor=(0.2, 0.05, 0.8, 0.8),
  bbox_transform=fig0.transFigure,
  loc="center"
)
axins.set_xlim(2.8, 3.0)
axins.set_ylim(0.0322, 0.0342)
# axins.yaxis.set_major_locator(MultipleLocator(0.0005))   # dy = 0.05
mark_inset(ax0, axins, loc1=2, loc2=4, fc="none", ec="0.45")
# axins.yaxis.tick_right()

if has_pde:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      axins.plot(pde_barycenter_data[i]["t"], pde_barycenter_data[i]["vx"], '-', lw=line_width, label=f"PDE-based ($f$={f:03d})", color = pde_colors[i])
  else:
    axins.plot(t_pde, y_pde, '-', lw=line_width, label="PDE-based")
if has_geo:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      axins.plot(geo_barycenter_data[i]["t"], geo_barycenter_data[i]["vx"], '-', lw=line_width, label=f"Geometric ($f$={f:03d})", color = geo_colors[i])
  else:
    axins.plot(t_geo, y_geo, '-', lw=line_width, label="Geometric")
if has_pro:
  if (args.frequency_sweep):
    for i,f in enumerate(f_values):
      axins.plot(pro_barycenter_data[i]["t"], pro_barycenter_data[i]["vx"], '-', lw=line_width, label=f"Projection ($f$={f:03d})", color = pro_colors[i])
  else:
    axins.plot(t_pro, y_pro, '-', lw=line_width, label="Projection")
# 
axins.plot(t_analytical, vx_analytical, ls=c1_ref_linestyles[0], c=c1_color, lw=line_width, alpha=c1_alpha, label="Analytic")




ax0.set_ylabel(r'Migration velocity [LT$^{-1}$]')
ax0.set_xlabel(r'Migration time [T]')
# ax0.legend(loc="upper left")
fig0.savefig(f'{args.save_figure_path}/migration_velocity.png',dpi=300)
fig0.savefig(f'{args.save_figure_path}/migration_velocity.svg',dpi=300)

ax0_proj.set_ylabel(r'Migration velocity [LT$^{-1}$]')
ax0_proj.set_xlabel(r'Migration time [T]')
# ax0_proj.legend(loc="upper left")
fig0_proj.savefig(f'{args.save_figure_path}/migration_velocity_proj.png',dpi=300)
fig0_proj.savefig(f'{args.save_figure_path}/migration_velocity_proj.svg',dpi=300)

plt.close()


################################
# Plot the sphericity evolution
################################

fig5 = plt.figure()
ax5 = fig5.add_subplot(111)

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)

# inset axes
# axins_geo = inset_axes(
#   ax5, width="40%", height="40%",
#   bbox_to_anchor=(0.8, -0.1, 0.8, 0.8),
#   bbox_transform=fig5.transFigure,
#   loc="center"
# )
# axins_geo.set_xlim(2.5, 3.0)
# axins_geo.set_ylim(0.78,0.79)
# mark_inset(ax5, axins_geo, loc1=2, loc2=3, fc="none", ec="0.45",zorder=3)
# axins_geo.yaxis.set_major_locator(MultipleLocator(0.005))
# axins_geo.yaxis.tick_right()

# axins_pde = inset_axes(
#   ax5, width="40%", height="40%",
#   bbox_to_anchor=(0.8, 0.3, 0.8, 0.8),
#   bbox_transform=fig5.transFigure,
#   loc="center"
# )
# axins_pde.set_xlim(2.5, 3.0)
# axins_pde.set_ylim(0.819, 0.822)
# mark_inset(ax5, axins_pde, loc1=2, loc2=3, fc="none", ec="0.45",zorder=3)
# axins_pde.yaxis.set_major_locator(MultipleLocator(0.005))
# axins_pde.yaxis.tick_right()

if has_pro:
  if (args.frequency_sweep):
    for i,freq in enumerate(f_values):
      df_mass_pro = pd.read_csv(pro_mass_geo_data[i], header=0, sep=r'\s+')
      # equivalent_surface_area = np.pi**(1/3) * (6*df_mass_pro["geometric_volume_fluid_1"])**(2/3)
      equivalent_surface_area = df_mass_pro['surface_fluid_1'][0]
      ax6.plot(df_mass_pro['time'], df_mass_pro['surface_fluid_1'],
               "-", label=f'Projection ($f$ = {freq:03d})', linewidth=line_width, color=pro_colors[i])

  else:
    equivalent_surface_area = np.pi**(1/3) * (6*df_mass_pro["geometric_volume_fluid_1"])**(2/3)
    ax6.plot(df_mass_pro['time'],
             df_mass_pro['surface_fluid_1'],
             "-", label=r'Projection', linewidth=line_width)
             
if has_pde:
  if (args.frequency_sweep):
    for i,freq in enumerate(f_values):
      df_mass_pde = pd.read_csv(pde_mass_geo_data[i], header=0, sep=r'\s+')
      # equivalent_surface_area = np.pi**(1/3) * (6*df_mass_pde["geometric_volume_fluid_1"])**(2/3)
      equivalent_surface_area = df_mass_pde['surface_fluid_1'][0]
      ax5.plot(df_mass_pde['time'], df_mass_pde['surface_fluid_1'],
               "-", label=f'PDE-based ($f$ = {freq:03d})', linewidth=line_width, color=pde_colors[i])
      # axins_pde.plot(df_mass_pde['time'], df_mass_pde['surface_fluid_1'],
              #  "-", label=f'PDE-based ($f$ = {freq:03d})', linewidth=line_width, color=pde_colors[i])
  else:
    equivalent_surface_area = np.pi**(1/3) * (6*df_mass_pde["geometric_volume_fluid_1"])**(2/3)
    ax5.plot(df_mass_pde['time'],
             df_mass_pde['surface_fluid_1'],
             "-", label=r'PDE-based', linewidth=line_width)
if has_geo:
  if (args.frequency_sweep):
    for i,freq in enumerate(f_values):
      df_mass_geo = pd.read_csv(geo_mass_geo_data[i], header=0, sep=r'\s+')
      # equivalent_surface_area = np.pi**(1/3) * (6*df_mass_geo["geometric_volume_fluid_1"])**(2/3)
      equivalent_surface_area = df_mass_geo['surface_fluid_1'][0]
      
      ax5.plot(df_mass_geo['time'], df_mass_geo['surface_fluid_1'],
               "-", label=f'Geometric ($f$ = {freq:03d})', linewidth=line_width, color=geo_colors[i])
      # axins_geo.plot(df_mass_geo['time'], df_mass_geo['surface_fluid_1'],
              #  "-", label=f'Geometric ($f$ = {freq:03d})', linewidth=line_width, color=geo_colors[i],zorder=1)
  else:
    equivalent_surface_area = np.pi**(1/3) * (6*df_mass_geo["geometric_volume_fluid_1"])**(2/3)
    ax5.plot(df_mass_geo['time'],
      df_mass_geo['surface_fluid_1'],
             "-", label=r'Geometric', linewidth=line_width)


ax5.set_ylabel(r'${A_\mathrm{droplet}} [L^2]$')
ax5.set_xlabel(r'Migration time [T]')
ax5.ticklabel_format(useOffset=False, style='plain', axis='y')
# ax5.set_ylim(0.95,1.01)
# ax5.legend(loc="upper left")
fig5.savefig(f'{args.save_figure_path}/migration-sphericity.png',dpi=300)
fig5.savefig(f'{args.save_figure_path}/migration-sphericity.svg',dpi=300)

# plt.show()
plt.close()

ax6.set_ylabel(r'${A_\mathrm{droplet}} [L^2]$')
ax6.set_xlabel(r'Migration time [T]')
ax6.ticklabel_format(useOffset=False, style='plain', axis='y')
# ax6.legend(loc="upper left")
fig6.savefig(f'{args.save_figure_path}/migration-sphericity-proj.png',dpi=300)
fig6.savefig(f'{args.save_figure_path}/migration-sphericity-proj.svg',dpi=300)

# plt.show()
plt.close()

#************************************************************
# Plot the contour of the bubble in the last frame
#************************************************************

# For slice extraction
origin = np.array([0,0,0])
normal = np.array([0,0,1])

fig_pro_radius = plt.figure()
ax_pro_radius = fig_pro_radius.add_subplot(111)

if has_pro:
  if (args.frequency_sweep):
    for i,freq in enumerate(f_values):
      if (freq == 1):
        continue
      x_pro_barycenter = pro_barycenter_data[i]["x"][-1]
      y_pro_barycenter = pro_barycenter_data[i]["y"][-1]
      z_pro_barycenter = pro_barycenter_data[i]["z"][-1]
      
      output_dir_pro = base_output_dir_pro.format(f=freq)

      # Get latest vtu file for slice extraction
      list_vtu = [f
                  for f in os.listdir(output_dir_pro)
                  if (f.endswith('.pvtu') and "interface" not in f and "sliced" not in f)]

      latest_file = os.path.join(output_dir_pro, max(list_vtu, key=extract_number))
      print(latest_file)

      # Slice
      slice_extractor.extract_slice(latest_file,origin,normal,output_dir_pro)
      output_files_list = os.listdir(output_dir_pro)
      sliced_list_vtu = [(output_dir_pro+"/"+x) for x in output_files_list if  ("sliced" in x)]
      latest_sliced_file = max(sliced_list_vtu, key=os.path.getctime)
      print("Opening file: ", latest_sliced_file)
      sim = pv.read(latest_sliced_file)
      sim.set_active_scalars("phase")
      contour_val = np.array([0.5])
      contours = sim.contour(contour_val)
      x_pro, y_pro, z_pro = contours.points[:, 0], contours.points[:, 1], contours.points[:, 2]

      # Store in dictionary
      pro_contour_data[i] = {
        "x": x_pro,
        "y": y_pro,
        "z": z_pro,
      }
      
      radius = np.empty(len(x_pro))
      angle = np.empty(len(x_pro))
      
      # Compute radius and angle
      for j in range(len(x_pro)):
          radius[j] = np.sqrt((x_pro_barycenter-x_pro[j])**2 + (y_pro_barycenter-y_pro[j])**2)
          angle[j] = np.arctan2((y_pro_barycenter-y_pro[j]),(x_pro_barycenter-x_pro[j]))
          angle[j] = np.where(angle[j] < 0, angle[j] + 2 * np.pi, angle[j])

      # Sort
      p = angle.argsort()
      
      # Set the x-axis limits
      ax_pro_radius.set_xlim(0, 2 * np.pi)

      # Define the tick locations (e.g., at 0, pi/2, pi, 3pi/2, 2pi)
      tick_locations = [0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
      ax_pro_radius.set_xticks(tick_locations)
  
      # Define the corresponding labels using LaTeX for pi
      tick_labels = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
      ax_pro_radius.set_xticklabels(tick_labels)

      ax_pro_radius.plot(angle[p], radius[p], '-', lw=line_width,label=f"Projection ($f$={f:03d})", color = pro_colors[i])
      
  else:
    # Get latest vtu file for slice extraction
    list_vtu = [f
                for f in os.listdir(output_dir_pro)
                if (f.endswith('.vtu') and "interface" not in f and "sliced" not in f)]
    # latest_file = os.path.join(output_dir_pro,max(list_vtu, key=extract_number))
    latest_file = os.path.join(output_dir_pro, max(list_vtu, key=extract_number))


    # Slice
    slice_extractor.extract_slice(latest_file,origin,normal,output_dir_pro)
    output_files_list = os.listdir(output_dir_pro)
    sliced_list_vtu = [(output_dir_pro+"/"+x) for x in output_files_list if  ("sliced" in x)]
    latest_sliced_file = max(sliced_list_vtu, key=os.path.getctime)
    print("Opening file: ", latest_sliced_file)
    sim = pv.read(latest_sliced_file)
    sim.set_active_scalars("phase")
    contour_val = np.array([0.5])
    contours = sim.contour(contour_val)
    x_pro, y_pro = contours.points[:, 0], contours.points[:, 1]

ax_pro_radius.plot(angle[p], 0.25*np.ones(len(angle[p])), '--k', lw=line_width)

ax_pro_radius.set_ylabel(r'Radius [L]')
ax_pro_radius.set_xlabel(r'Angle [rad]')
ax_pro_radius.set_ylim(0.24, 0.26)

fig_pro_radius.savefig(f'{args.save_figure_path}/radius_pro.png',dpi=300)
fig_pro_radius.savefig(f'{args.save_figure_path}/radius_pro.svg',dpi=300)
plt.close()

fig_geo_radius = plt.figure()
ax_geo_radius = fig_geo_radius.add_subplot(111)

if has_geo:
  if (args.frequency_sweep):
    for i,freq in enumerate(f_values):
        
      x_geo_barycenter = geo_barycenter_data[i]["x"][-1]
      y_geo_barycenter = geo_barycenter_data[i]["y"][-1]
      z_geo_barycenter = geo_barycenter_data[i]["z"][-1]
      
      output_dir_geo = base_output_dir_geo.format(f=freq)

      # Get latest vtu file for slice extraction
      list_vtu = [f
                  for f in os.listdir(output_dir_geo)
                  if (f.endswith('.pvtu') and "interface" not in f and "sliced" not in f)]
      # latest_file = os.path.join(output_dir_geo,max(list_vtu, key=extract_number))
      # latest_file = os.path.join(output_dir_geo,list_vtu[-1])
      latest_file = os.path.join(output_dir_geo, max(list_vtu, key=extract_number))


      # Slice
      slice_extractor.extract_slice(latest_file,origin,normal,output_dir_geo)
      output_files_list = os.listdir(output_dir_geo)
      sliced_list_vtu = [(output_dir_geo+"/"+x) for x in output_files_list if  ("sliced" in x)]
      latest_sliced_file = max(sliced_list_vtu, key=os.path.getctime)
      print("Opening file: ", latest_sliced_file)
      sim = pv.read(latest_sliced_file)
      sim.set_active_scalars("phase")
      contour_val = np.array([0.5])
      contours = sim.contour(contour_val)
      x_geo, y_geo, z_geo = contours.points[:, 0], contours.points[:, 1], contours.points[:, 2]

      # Store in dictionary
      geo_contour_data[i] = {
        "x": x_geo,
        "y": y_geo,
        "z": z_geo,      
      }
      
      radius = np.empty(len(x_geo))
      angle = np.empty(len(x_geo))
      
      # Compute radius and angle
      for j in range(len(x_geo)):
          radius[j] = np.sqrt((x_geo_barycenter-x_geo[j])**2 + (y_geo_barycenter-y_geo[j])**2)
          angle[j] = np.arctan2((y_geo_barycenter-y_geo[j]),(x_geo_barycenter-x_geo[j]))
          angle[j] = np.where(angle[j] < 0, angle[j] + 2 * np.pi, angle[j])

      # Sort
      p = angle.argsort()
      
      # Set the x-axis limits
      ax_geo_radius.set_xlim(0, 2 * np.pi)

      # Define the tick locations (e.g., at 0, pi/2, pi, 3pi/2, 2pi)
      tick_locations = [0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
      ax_geo_radius.set_xticks(tick_locations)
  
      # Define the corresponding labels using LaTeX for pi
      tick_labels = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
      ax_geo_radius.set_xticklabels(tick_labels)

      ax_geo_radius.plot(angle[p], radius[p], '-', lw=line_width,label=f"Geometric ($f$={f:03d})", color = geo_colors[i])
      

  else:
    # Get latest vtu file for slice extraction
    list_vtu = [f
                for f in os.listdir(output_dir_geo)
                if (f.endswith('.vtu') and "interface" not in f and "sliced" not in f)]
    latest_file = os.path.join(output_dir_geo,max(list_vtu, key=extract_number))

    # Slice
    slice_extractor.extract_slice(latest_file,origin,normal,output_dir_geo)
    output_files_list = os.listdir(output_dir_geo)
    sliced_list_vtu = [(output_dir_geo+"/"+x) for x in output_files_list if  ("sliced" in x)]
    latest_sliced_file = max(sliced_list_vtu, key=os.path.getctime)
    print("Opening file: ", latest_sliced_file)
    sim = pv.read(latest_sliced_file)
    sim.set_active_scalars("phase")
    contour_val = np.array([0.5])
    contours = sim.contour(contour_val)
    x_geo, y_geo = contours.points[:, 0], contours.points[:, 1]

ax_geo_radius.plot(angle[p],0.25*np.ones(len(angle[p])), '--k', lw=line_width)

ax_geo_radius.set_ylabel(r'Radius [L]')
ax_geo_radius.set_xlabel(r'Angle [rad]')
ax_geo_radius.set_ylim(0.24, 0.26)

fig_geo_radius.savefig(f'{args.save_figure_path}/radius_geo.png',dpi=300)
fig_geo_radius.savefig(f'{args.save_figure_path}/radius_geo.svg',dpi=300)
plt.close()

fig_pde_radius = plt.figure()
ax_pde_radius = fig_pde_radius.add_subplot(111)

if has_pde:
  if (args.frequency_sweep):
    for i,freq in enumerate(f_values):
        
      x_pde_barycenter = pde_barycenter_data[i]["x"][-1]
      y_pde_barycenter = pde_barycenter_data[i]["y"][-1]
      z_pde_barycenter = pde_barycenter_data[i]["z"][-1]
      
      output_dir_pde = base_output_dir_pde.format(f=freq)

      # Get latest vtu file for slice extraction
      list_vtu = [f
                  for f in os.listdir(output_dir_pde)
                  if (f.endswith('.pvtu') and "interface" not in f and "sliced" not in f)]
      # latest_file = os.path.join(output_dir_pde,max(list_vtu, key=extract_number))
      # latest_file = os.path.join(output_dir_pde,list_vtu[-1])
      latest_file = os.path.join(output_dir_pde, max(list_vtu, key=extract_number))


      # Slice
      slice_extractor.extract_slice(latest_file,origin,normal,output_dir_pde)
      output_files_list = os.listdir(output_dir_pde)
      sliced_list_vtu = [(output_dir_pde+"/"+x) for x in output_files_list if  ("sliced" in x)]
      latest_sliced_file = max(sliced_list_vtu, key=os.path.getctime)
      print("Opening file: ", latest_sliced_file)
      sim = pv.read(latest_sliced_file)
      sim.set_active_scalars("phase")
      contour_val = np.array([0.5])
      contours = sim.contour(contour_val)
      x_pde, y_pde, z_pde = contours.points[:, 0], contours.points[:, 1], contours.points[:, 2]

      # Store in dictionary
      pde_contour_data[i] = {
        "x": x_pde,
        "y": y_pde,
        "z": z_pde,
      }
      
      radius = np.empty(len(x_pde))
      angle = np.empty(len(x_pde))
      
      # Compute radius and angle
      for j in range(len(x_pde)):
          radius[j] = np.sqrt((x_pde_barycenter-x_pde[j])**2 + (y_pde_barycenter-y_pde[j])**2)
          angle[j] = np.arctan2((y_pde_barycenter-y_pde[j]),(x_pde_barycenter-x_pde[j]))
          angle[j] = np.where(angle[j] < 0, angle[j] + 2 * np.pi, angle[j])

      # Sort
      p = angle.argsort()
      
      # Set the x-axis limits
      ax_pde_radius.set_xlim(0, 2 * np.pi)

      # Define the tick locations (e.g., at 0, pi/2, pi, 3pi/2, 2pi)
      tick_locations = [0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi]
      ax_pde_radius.set_xticks(tick_locations)
  
      # Define the corresponding labels using LaTeX for pi
      tick_labels = [r'$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$']
      ax_pde_radius.set_xticklabels(tick_labels)

      ax_pde_radius.plot(angle[p], radius[p], '-', lw=line_width,label=f"Projection ($f$={f:03d})", color = pde_colors[i])
  else:
    # Get latest vtu file for slice extraction
    list_vtu = [f
                for f in os.listdir(output_dir_pde)
                if (f.endswith('.vtu') and "interface" not in f and "sliced" not in f)]
    latest_file = os.path.join(output_dir_pde,max(list_vtu, key=extract_number))

    # Slice
    slice_extractor.extract_slice(latest_file,origin,normal,output_dir_pde)
    output_files_list = os.listdir(output_dir_pde)
    sliced_list_vtu = [(output_dir_pde+"/"+x) for x in output_files_list if  ("sliced" in x)]
    latest_sliced_file = max(sliced_list_vtu, key=os.path.getctime)
    print("Opening file: ", latest_sliced_file)
    sim = pv.read(latest_sliced_file)
    sim.set_active_scalars("phase")
    contour_val = np.array([0.5])
    contours = sim.contour(contour_val)
    x_pde, y_pde = contours.points[:, 0], contours.points[:, 1]

ax_pde_radius.plot(angle[p], 0.25*np.ones(len(angle[p])), '--k', lw=line_width)

ax_pde_radius.set_ylabel(r'Radius [L]')
ax_pde_radius.set_xlabel(r'Angle [rad]')
ax_pde_radius.set_ylim(0.24, 0.26)

fig_pde_radius.savefig(f'{args.save_figure_path}/radius_pde.png',dpi=300)
fig_pde_radius.savefig(f'{args.save_figure_path}/radius_pde.svg',dpi=300)
plt.close()

# Get path to CFL evolution file
base_dir_pde = os.path.dirname(args.pde.rstrip("/")) + "/"
base_dir_geo = os.path.dirname(args.geo.rstrip("/")) + "/"
base_dir_pro = os.path.dirname(args.pro.rstrip("/")) + "/"
base_dir_pde = base_format(base_dir_pde)
base_dir_geo = base_format(base_dir_geo)
base_dir_pro = base_format(base_dir_pro)

# Load CFL evolution data for each method
cfl_pde = load_cfl_evolution(has_pde, base_dir_pde, f_values)
cfl_geo = load_cfl_evolution(has_geo, base_dir_geo, f_values)
cfl_pro = load_cfl_evolution(has_pro, base_dir_pro, f_values)

# Plot CFL evolution
plot_cfl_evolution(cfl_pde, cfl_geo, cfl_pro, args.save_figure_path, f_values, line_width, pde_colors, geo_colors, pro_colors)

# Keep only the defined values
# x, y, label, fig_name, color_list = [], [], [], [], []

# if has_pde:
#   if (args.frequency_sweep):
#     for i,freq in enumerate(f_values):
#       x.append(pde_contour_data[i]['x'])
#       y.append(pde_contour_data[i]['y'])
#       label.append(f'PDE-based ($f$ = {freq:03d})')
#       fig_name.append(f'pde_f{freq:03d}')
#       color_list.append(pde_colors[i])
#   else:
#     x.append(x_pde)
#     y.append(y_pde)
#     label.append('PDE-based')
#     fig_name.append('pde')
# if has_geo:
#   if (args.frequency_sweep):
#     for i,freq in enumerate(f_values):
#       x.append(geo_contour_data[i]['x'])
#       y.append(geo_contour_data[i]['y'])
#       label.append(f'Geometric ($f$ = {freq:03d})')
#       fig_name.append(f'geo_f{freq:03d}')
#       color_list.append(geo_colors[i])
#   else:
#     x.append(x_geo)
#     y.append(y_geo)
#     label.append('Geometric')
#     fig_name.append('geo')
# if has_pro:
#   if (args.frequency_sweep):
#     for i,freq in enumerate(f_values):
#       if (freq == 1):
#         continue
#       x.append(pro_contour_data[i]['x'])
#       y.append(pro_contour_data[i]['y'])
#       label.append(f'Projection ($f$ = {freq:03d})')
#       fig_name.append(f'pro_f{freq:03d}')
#       color_list.append(pro_colors[i])
#   else:
#     x.append(x_pro)
#     y.append(y_pro)
#     label.append('Projection')
#     fig_name.append('pro')

# # Plot bubble contour comparison between the three methods (only if more than 1 method)
# if n > 1:
#   str_method = ["pde", "geo", "pro"]
#   method_id = 0
#   i_start = 0
#   while (i_start <= n-len(f_values)):
#     fig2 = plt.figure()
#     ax2= fig2.add_subplot(111)

#     for i in range(i_start,i_start+len(f_values)):
#       ax2.scatter(x[i], y[i], s=4, marker="s", color=color_list[i], label=label[i], linewidths=line_width)

#     ax2.set_xlabel(r'x [L]')
#     ax2.set_ylabel(r'y [L]')
#     # ax2.legend(markerscale=2, scatterpoints=20)
#     ax2.grid( which='major', color='grey', linestyle='--')
#     # ax2.set_xlim([0.2,0.8])
#     # ax2.set_ylim([1.15,1.75])
#     i_start += len(f_values)
#     fig2.savefig(f"{args.save_figure_path}/bubble-contour-" + str_method[method_id] + ".png",dpi=300)
#     method_id += 1
#     plt.close()
#       # plt.show()

