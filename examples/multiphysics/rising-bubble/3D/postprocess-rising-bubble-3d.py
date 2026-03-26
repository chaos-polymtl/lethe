# SPDX-FileCopyrightText: Copyright (c) 2026 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

#############################################################################
"""
Postprocessing code for rising-bubble example
"""
#############################################################################

'''Importing Libraries'''
from rising_bubble_utils import *
import argparse
import os
#############################################################################
def parse_args():
  '''
  Parse arguments of the script
  :return: Parsed arguments
  '''
  parser = argparse.ArgumentParser(description='Arguments for the validation of the 3D rising bubble benchmark')
  parser.add_argument("-p", "--pro", type=str, help="Path to the output folder for the projecting results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=False)
  parser.add_argument("-g", "--geo", type=str, help="Path to the output folder for the geometric results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=False)
  parser.add_argument("-a", "--pde", type=str, help="Path to the output folder for the PDE-based results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=False)
  parser.add_argument("-n", "--none", type=str, help="Path to the output folder for the none results. This is the folder that contains the results of the simulation (.vtu, .pvtu, .dat and .pvd files)", required=False)
  parser.add_argument("-l", "--lethe_path", type=str, help="Path to the Lethe's main folder", required=True)
  parser.add_argument("-sf", "--save_figure_path", type=str, help="Path to where figure should be saved", required=False, default="./")
  parser.add_argument("-cfl", "--cfl_evolution", type=bool, help="Set to true to plot CFL evolution figures", required=False, default=False)
  return parser.parse_args()

#******************************
# Global style and constants
#******************************
pro_colors = ['#006d2c','#31a354','#74c476','#bae4b3']
pde_colors = ['#54278f','#756bb1','#9e9ac8','#cbc9e2']
geo_colors = ['#a63603','#e6550d','#fd8d3c','#fdbe85']
none_colors = ['#808080']

c1_ref_linestyles = ['--',(0, (3, 1, 1, 1, 1, 1)), '-.', ':']
c1_ref_markers = ['', '', '', '']
c1_ref_label = ["DROPS (2014)", "FeatFlow (2019)", "NaSt3DGPF (2014)", "OpenFOAM (2014)"]
c1_alpha = 1.0
c1_color = 'k'

lethe_linestyle = "-"

line_width = 3.0

case_number = 1         # Only case #1 is post-processed here
f_values = [1, 10, 100] # Frequencies studied

""" MAIN """
def main():
  setup_plot_params()
  args = parse_args()

  # Check post-processed methods
  has_pde = args.pde is not None
  has_geo = args.geo is not None
  has_pro = args.pro is not None
  has_none = args.none is not None


  # Create save folder if it does not exist
  os.makedirs(args.save_figure_path, exist_ok=True)

  # Import VTU slicer
  slice_extractor = import_slice_extractor(args.lethe_path)

  # Get base output directory names for each method
  base_output_dir_pde = base_format(args.pde)
  base_output_dir_geo = base_format(args.geo)
  base_output_dir_pro = base_format(args.pro)
  base_output_dir_none = base_format(args.none)


  # Load Lethe results
  pde_bary, pde_mass = load_method_data(has_pde, base_output_dir_pde, f_values)
  geo_bary, geo_mass = load_method_data(has_geo, base_output_dir_geo, f_values)
  pro_bary, pro_mass = load_method_data(has_pro, base_output_dir_pro, f_values)
  none_bary, none_mass = load_method_data(has_none, base_output_dir_none, [0])


  # Load reference results
  refs = read_refs()

  # Generate legends for figures
  make_legend_no_ref(f'{args.save_figure_path}/legend_no_ref_c{case_number}', f_values,  c1_ref_linestyles,
                          c1_alpha, lethe_linestyle, pde_colors, geo_colors, pro_colors)  
  make_legend_no_featflow(f'{args.save_figure_path}/legend_no_featflow_c{case_number}', f_values,  c1_ref_linestyles,
                          c1_alpha, lethe_linestyle, pde_colors, geo_colors, pro_colors)
  make_legend_with_featflow(f'{args.save_figure_path}/legend_with_featflow_c{case_number}', f_values, c1_ref_linestyles,
                            c1_ref_markers, c1_alpha, lethe_linestyle, pde_colors, geo_colors, pro_colors)


  # Generate plots for analyses
  plot_barycenter(pde_bary, geo_bary, pro_bary, none_bary, refs, args.save_figure_path, f_values, case_number, line_width,
                  c1_ref_linestyles, c1_color, c1_alpha, pde_colors, geo_colors, pro_colors, none_colors)
  plot_rise_velocity(pde_bary, geo_bary, pro_bary, none_bary, refs, args.save_figure_path, f_values, case_number, line_width,
                     c1_ref_linestyles, c1_color, c1_alpha, pde_colors, geo_colors, pro_colors, none_colors)
  pde_c, geo_c, pro_c = extract_contours(slice_extractor, has_pde, has_geo, has_pro,
                                         base_output_dir_pde, base_output_dir_geo, base_output_dir_pro, f_values)
  plot_global_mass_conservation(pde_mass, geo_mass, pro_mass, args.save_figure_path, f_values, case_number, line_width,
                                pde_colors, geo_colors, pro_colors)
  plot_geometric_mass_conservation(pde_mass, geo_mass, pro_mass, none_mass, args.save_figure_path, f_values, case_number, line_width,
                                   pde_colors, geo_colors, pro_colors, none_colors)
  plot_sphericity(pde_mass, geo_mass, pro_mass, none_mass, refs, args.save_figure_path, f_values, case_number, line_width,
                  c1_ref_linestyles, c1_color, c1_alpha, pde_colors, geo_colors, pro_colors, none_colors)
  plot_contour_sets(pde_c, geo_c, pro_c, args.save_figure_path, f_values, case_number, line_width)
  

  if (args.cfl_evolution):
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
    plot_cfl_evolution(cfl_pde, cfl_geo, cfl_pro, args.save_figure_path, f_values, case_number, line_width, pde_colors, geo_colors, pro_colors)

if __name__ == "__main__":
  main()