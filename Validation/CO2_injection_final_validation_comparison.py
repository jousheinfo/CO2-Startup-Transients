# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 10:57:49 2025

@author: jpauya1
"""


import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from t2listing import *
import CoolProp.CoolProp as CP

def save_figure_png_and_eps(fig, filename, dpi=300):
    """Save a Matplotlib figure in both PNG and EPS formats."""
    png_path = Path(filename)
    eps_path = png_path.with_suffix(".eps")

    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(eps_path, format="eps", dpi=dpi, bbox_inches="tight")

    print(f"Saved: {png_path}")
    print(f"Saved: {eps_path}")
    
#Saturation Line
data_SL = pd.read_csv('CO2_Saturation_Line_NIST.txt',delim_whitespace=True,header=None,skiprows=(3))
data_SL.columns = ['Temperature','Pressure']
print(data_SL)

#Setting the x- and y-axis for SL
T_SL = data_SL['Temperature']
P_SL = data_SL['Pressure']

# ---------------------------------------------------------------------------

#T2WELL

#General output file
lst = t2listing('CO2_injection_no_caprock_ECO2N_final_test_no_salinity_pytough.listing')

#----------------------------------------------------------------------------

#Wellhead
wellhead_name = '*wa 1'
well_head = wellhead_name

#Bottomhole
above_bottomhole_name = '  a80'
above_bottom_hole = above_bottomhole_name

#Bottomhole
bottomhole_name_1 = '  a81'
bottom_hole_1 = bottomhole_name_1

bottomhole_name_2 = '  a82'
bottom_hole_2 = bottomhole_name_2

bottomhole_name_3 = '  a83'
bottom_hole_3 = bottomhole_name_3

bottomhole_name_4 = '  a84'
bottom_hole_4 = bottomhole_name_4

bottomhole_name_5 = '  a85'
bottom_hole_5 = bottomhole_name_5

bottomhole_name_6 = '  a86'
bottom_hole_6 = bottomhole_name_6

bottomhole_name_7 = '  a87'
bottom_hole_7 = bottomhole_name_7

bottomhole_name_8 = '  a88'
bottom_hole_8 = bottomhole_name_8

bottomhole_name_9 = '  a89'
bottom_hole_9 = bottomhole_name_9

bottomhole_name_10 = '  a90'
bottom_hole_10 = bottomhole_name_10

# #Reservoir
# reservoir_name_1_1 = '  b59'
# reservoir_1_1 = reservoir_name_1_1

# reservoir_name_10_1 = '  b68'
# reservoir_10_1 = reservoir_name_10_1

#----------------------------------------------------------------------------

#Defining the pressure and temperature values

#Wellhead
[(twh,pwh), (twh, tempwh), (twh, dgwh)] = lst.history([('e', well_head, 'P'), ('e', well_head, 'T'), ('e', well_head, 'DG')]) #e stands for 'element'

# Above bottomhole
[(tabh, pabh), (tabh, tempabh), (tabh,dgabh)] = lst.history([('e', above_bottom_hole, 'P'), ('e', above_bottom_hole, 'T'), ('e', above_bottom_hole, 'DG')]) #e stands for 'element'

#Bottomhole
[(tbh_1, pbh_1), (tbh_1, tempbh_1), (tbh_1,dgbh_1)] = lst.history([('e', bottom_hole_1, 'P'), ('e', bottom_hole_1, 'T'), ('e', bottom_hole_1, 'DG')]) #e stands for 'element'
[(tbh_2, pbh_2), (tbh_2, tempbh_2), (tbh_2,dgbh_2)] = lst.history([('e', bottom_hole_2, 'P'), ('e', bottom_hole_2, 'T'), ('e', bottom_hole_2, 'DG')]) #e stands for 'element'
[(tbh_3, pbh_3), (tbh_3, tempbh_3), (tbh_3,dgbh_3)] = lst.history([('e', bottom_hole_3, 'P'), ('e', bottom_hole_3, 'T'), ('e', bottom_hole_3, 'DG')]) #e stands for 'element'
[(tbh_4, pbh_4), (tbh_4, tempbh_4), (tbh_4,dgbh_4)] = lst.history([('e', bottom_hole_4, 'P'), ('e', bottom_hole_4, 'T'), ('e', bottom_hole_4, 'DG')]) #e stands for 'element'
[(tbh_5, pbh_5), (tbh_5, tempbh_5), (tbh_5,dgbh_5)] = lst.history([('e', bottom_hole_5, 'P'), ('e', bottom_hole_5, 'T'), ('e', bottom_hole_5, 'DG')]) #e stands for 'element'
[(tbh_6, pbh_6), (tbh_6, tempbh_6), (tbh_6,dgbh_6)] = lst.history([('e', bottom_hole_6, 'P'), ('e', bottom_hole_6, 'T'), ('e', bottom_hole_6, 'DG')]) #e stands for 'element'
[(tbh_7, pbh_7), (tbh_7, tempbh_7), (tbh_7,dgbh_7)] = lst.history([('e', bottom_hole_7, 'P'), ('e', bottom_hole_7, 'T'), ('e', bottom_hole_7, 'DG')]) #e stands for 'element'
[(tbh_8, pbh_8), (tbh_8, tempbh_8), (tbh_8,dgbh_8)] = lst.history([('e', bottom_hole_8, 'P'), ('e', bottom_hole_8, 'T'), ('e', bottom_hole_8, 'DG')]) #e stands for 'element'
[(tbh_9, pbh_9), (tbh_9, tempbh_9), (tbh_9,dgbh_9)] = lst.history([('e', bottom_hole_9, 'P'), ('e', bottom_hole_9, 'T'), ('e', bottom_hole_9, 'DG')]) #e stands for 'element'
[(tbh_10, pbh_10), (tbh_10, tempbh_10), (tbh_10,dgbh_10)] = lst.history([('e', bottom_hole_10, 'P'), ('e', bottom_hole_10, 'T'), ('e', bottom_hole_10, 'DG')]) #e stands for 'element'

tbh = (tbh_1 + tbh_2 + tbh_3 + tbh_4 + tbh_5 + tbh_6 + tbh_7 + tbh_8 + tbh_9 + tbh_10) / 10
pbh = (pbh_1 + pbh_2 + pbh_3 + pbh_4 + pbh_5 + pbh_6 + pbh_7 + pbh_8 + pbh_9 + pbh_10) / 10
tempbh = (tempbh_1 + tempbh_2 + tempbh_3 + tempbh_4 + tempbh_5 + tempbh_6 + tempbh_7 + tempbh_8 + tempbh_9 + tempbh_10) / 10
dgbh = (dgbh_1 + dgbh_2 + dgbh_3 + dgbh_4 + dgbh_5 + dgbh_6 + dgbh_7 + dgbh_8 + dgbh_9 + dgbh_10) / 10

#Reservoir
# [(tres_1_1, pres_1_1), (tres_1_1, tempres_1_1), (tres,dgres_1_1)] = lst.history([('e', reservoir_1_1, 'P'), ('e', reservoir_1_1, 'T'), ('e', reservoir_1_1, 'DG')]) #e stands for 'element'
# [(tres_10_1, pres_10_1), (tres_10_1, tempres_10_1), (tres,dgres_10_1)] = lst.history([('e', reservoir_10_1, 'P'), ('e', reservoir_10_1, 'T'), ('e', reservoir_10_1, 'DG')]) #e stands for 'element'

#----------------------------------------------------------------------------

#Saturation Line
data_SL = pd.read_csv('CO2_Saturation_Line_NIST.txt', delim_whitespace=True,header=None,skiprows=(3))
data_SL.columns = ['Temperature','Pressure']
#print(data_SL)

#Setting the x- and y-axis for SL
T_SL = data_SL['Temperature']
P_SL = data_SL['Pressure']

#----------------------------------------------------------------------------

fflow_data = pd.read_csv('FFlow',sep=',',header=None,skiprows=(3))
fflow_data.columns = ['Time','Depth','Fliq','Fgas','Vliq','Vgas','Umix']
above_perforations = fflow_data[fflow_data['Depth'] == 0.1470003019677E+04]

x = above_perforations['Time']
x_above_perforations = x
y_above_perforations_umix = above_perforations['Umix']

#----------------------------------------------------------------------------

#Plotting pressure vs time

#Steady State reached after 6 hours

# fig, axs = plt.subplots(1, 3, figsize=(17,6))
# #ancho, alto

# axs[0].plot(twh[:42] / 3600, pwh[:42] * 0.00001, color='blue', label='wellhead')
# axs[0].plot(tbh[:42] / 3600, pbh[:42] * 0.00001, color='red', label='bottomhole')
# # axs[0].plot(tres_1_1[:42] / 3600, pres_1_1[:42] * 0.00001, color='purple', label='1st reservoir')
# # axs[0].plot(tres_10_1[:42] / 3600, pres_10_1[:42] * 0.00001, color='red', label='1st reservoir')
# axs[0].axhline(y=P_SL.iloc[-1], xmin=0, xmax=twh[42]/3600, linewidth=2, color='black', linestyle='--')
# axs[0].legend(bbox_to_anchor=(0, 1.06, 3.2, 1.7), loc="lower center", borderaxespad=0, ncol=3)
# axs[0].set(xlabel ='Time (hours)', ylabel = 'Pressure (bar)')
# #axs[0].set_yticks([24, 56, 88, 120, 152, 184, 216, 248, 280])

# #Arrow to indicate
# axs[0].annotate(
# # Label and coordinate
# 'CO2 Pcrit = 73.773 bar', xy=(1.5, 74), xytext=(1.5, 100),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

# #Arrow to indicate
# axs[0].annotate(
# # Label and coordinate
# 'Initial WHP = 45 bar', xy=(0.05, 45), xytext=(1.5, 45),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

# #Arrow to indicate
# axs[0].annotate(
# # Label and coordinate
# 'Initial Pres = 243.8075 bar', xy=(0.05, 243.8075), xytext=(1.3, 220),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

# axs[1].plot(twh[:42] / 3600, tempwh[:42], color='blue', label='wellhead')
# axs[1].plot(tbh[:42] / 3600, tempbh[:42], color='red', label='bottomhole')
# # axs[1].plot(tres_1_1[:42] / 3600, tempres_1_1[:42], color='purple', label='1st reservoir')
# # axs[1].plot(tres_10_1[:42] / 3600, tempres_10_1[:42], color='red', label='1st reservoir')
# axs[1].axhline(y=T_SL.iloc[-1], xmin=0, xmax=twh[42]/3600, linewidth=2, color='black', linestyle='--')
# axs[1].set(xlabel ='Time (hours)', ylabel = 'Temperature (°C)')
# #axs[1].set_yticks([27, 32, 37, 42, 47, 52, 57, 62])

# #Arrow to indicate
# axs[1].annotate(
# # Label and coordinate
# 'CO2 Tcrit = 30.9782 °C', xy=(0.5, 31.3), xytext=(1.5, 33),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

# #Arrow to indicate
# axs[1].annotate(
# # Label and coordinate
# 'WHT = 35 °C', xy=(0.05, 35.3), xytext=(1.2, 38),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

# #Arrow to indicate
# axs[1].annotate(
# # Label and coordinate
# 'Initial Tres = 60 °C', xy=(0.05, 60), xytext=(1.4, 57),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

# axs[2].plot(twh[:42] / 3600, dgwh[:42], color='blue', label='wellhead')
# axs[2].plot(tbh[:42] / 3600, dgbh[:42], color='red', label='bottomhole')
# # axs[2].plot(tres_1_1[:42] / 3600, dgres_1_1[:42], color='purple', label='1st reservoir')
# # axs[2].plot(tres_10_1[:42] / 3600, dgres_10_1[:42], color='red', label='1st reservoir')
# axs[2].set(xlabel ='Time (hours)', ylabel = 'Gas density (kg/m3)')
# #axs[2].set_yticks([36, 101, 166, 231, 296, 361, 426, 491, 556, 621, 686, 751, 816, 881])

# plt.show()

# ---------------------------------------------------------------------------

# CO2LINK

# Reading the CSV file

n_cols = pd.read_csv('WH_Density_Plot_CO2LINK_no_caprock_no_salinity_final.csv', nrows=1, header=None).shape[1]
df_WH_Density = pd.read_csv('WH_Density_Plot_CO2LINK_no_caprock_no_salinity_final.csv', header=None, skiprows=2, usecols=[0, n_cols-2])
df_WH_Density.columns = ['Time (s)', 'Density_Ledaflow (kg/m3)' ]
time_WH_Density = df_WH_Density['Time (s)']
WH_Density_Ledaflow = df_WH_Density['Density_Ledaflow (kg/m3)']

df_WHP = pd.read_csv('Well_WHP_Plot_CO2LINK_no_caprock_no_salinity_final.csv', header=None, skiprows=2, usecols=[0, 1])
df_WHP.columns = ['Time (s)', 'WHP_Ledaflow (bar)' ]
time_WHP = df_WHP['Time (s)']
WHP_Ledaflow = df_WHP['WHP_Ledaflow (bar)']

# Calculate wellhead density using CoolProp for pure CO2
# Constant wellhead temperature: 35°C = 308.15 K
T_wellhead_K = 35.0 + 273.15  # Convert to Kelvin

# Calculate density for each pressure point
WH_Density_CoolProp = []
for pressure_bar in WHP_Ledaflow:
    pressure_Pa = pressure_bar * 1e5  # Convert bar to Pa
    try:
        # Calculate density using CoolProp for pure CO2
        density = CP.PropsSI('D', 'T', T_wellhead_K, 'P', pressure_Pa, 'CO2')
        WH_Density_CoolProp.append(density)
    except:
        # Handle any CoolProp errors
        WH_Density_CoolProp.append(float('nan'))

WH_Density_CoolProp = pd.Series(WH_Density_CoolProp)

# df_BH_Density = pd.read_csv('BH_Density_Plot_CO2LINK_no_caprock_no_salinity_final.csv', header=None, skiprows=2, usecols=[0, n_cols-2])
# df_BH_Density.columns = ['Time (s)', 'Density_Ledaflow (bar)' ]
# time_BH_Density = df_BH_Density['Time (s)']
# BH_Density_Ledaflow = df_BH_Density['Density_Ledaflow (bar)']

df_BHP = pd.read_csv('Well_BHP_Plot_CO2LINK_no_caprock_no_salinity_final.csv', header=None, skiprows=2)
df_BHP.columns = ['Time (s)', 'BHP_GEM (kPa)', 'BHP_Ledaflow (kPa)' ]
time_BHP = df_BHP['Time (s)']
BHP_Ledaflow = df_BHP['BHP_Ledaflow (kPa)']
BHP_GEM = df_BHP['BHP_GEM (kPa)']

df_BHT = pd.read_csv('Well_BHT_Plot_CO2LINK_no_caprock_no_salinity_final.csv', header=None, skiprows=2)
df_BHT.columns = ['Time (s)', 'BHT_GEM (C)', 'BHT_Ledaflow (C)']
time_BHT = df_BHT['Time (s)']
BHT_Ledaflow = df_BHT['BHT_Ledaflow (C)']
BHT_GEM = df_BHT['BHT_GEM (C)']

df_BHT = pd.read_csv('Well_BHT_Plot_CO2LINK_no_caprock_no_salinity_final.csv', header=None, skiprows=2)
df_BHT.columns = ['Time (s)', 'BHT_GEM (C)', 'BHT_Ledaflow (C)']
time_BHT = df_BHT['Time (s)']
BHT_Ledaflow = df_BHT['BHT_Ledaflow (C)']
BHT_GEM = df_BHT['BHT_GEM (C)']

df_above_perforations_velocity = pd.read_csv('Well_above_perforations_velocity_Plot_CO2LINK_no_caprock_no_salinity_final.csv', header=None, skiprows=2, usecols=[0, 1])
df_above_perforations_velocity.columns = ['Time (s)', 'above_perforations_gas_velocity_Ledaflow (m/s)']
time_above_perforations_velocity = df_above_perforations_velocity['Time (s)']
above_perforations_velocity_Ledaflow = df_above_perforations_velocity['above_perforations_gas_velocity_Ledaflow (m/s)']

# df_BHP_Tres_65C_Tres_20C = pd.read_csv('Well_Bottom-hole_Pressure_PlotData_Tres-65C_Tinj-20C.csv', header=None, skiprows=2)
# df_BHP_Tres_65C_Tres_20C.columns = ['Time (s)', 'BHP_GEM (kPa)', 'BHP_Ledaflow (kPa)' ]
# time_BHP_Tres_65C_Tres_20C = df_BHP_Tres_65C_Tres_20C['Time (s)']
# BHP_Tres_65C_Tres_20C = df_BHP_Tres_65C_Tres_20C['BHP_Ledaflow (kPa)']

# df_BHT_Tres_65C_Tres_20C = pd.read_csv('Well_Bottom-hole_Temperature_PlotData_Tres-65C_Tinj-20C.csv', header=None, skiprows=2)
# df_BHT_Tres_65C_Tres_20C.columns = ['Time (s)', 'BHT_GEM (C)', 'BHT_Ledaflow (C)']
# time_BHT_Tres_65C_Tres_20C = df_BHT_Tres_65C_Tres_20C['Time (s)']
# BHT_Tres_65C_Tres_20C = df_BHT_Tres_65C_Tres_20C['BHT_Ledaflow (C)']

# df_BHP_Tres_55C_Tres_20C = pd.read_csv('Well_Bottom-hole_Pressure_PlotData_Tres-55C_Tinj-20C.csv', header=None, skiprows=2)
# df_BHP_Tres_55C_Tres_20C.columns = ['Time (s)', 'BHP_GEM (kPa)', 'BHP_Ledaflow (kPa)' ]
# time_BHP_Tres_55C_Tres_20C = df_BHP_Tres_55C_Tres_20C['Time (s)']
# BHP_Tres_55C_Tres_20C = df_BHP_Tres_55C_Tres_20C['BHP_Ledaflow (kPa)']

# df_BHT_Tres_55C_Tres_20C = pd.read_csv('Well_Bottom-hole_Temperature_PlotData_Tres-55C_Tinj-20C.csv', header=None, skiprows=2)
# df_BHT_Tres_55C_Tres_20C.columns = ['Time (s)', 'BHT_GEM (C)', 'BHT_Ledaflow (C)']
# time_BHT_Tres_55C_Tres_20C = df_BHT_Tres_55C_Tres_20C['Time (s)']
# BHT_Tres_55C_Tres_20C = df_BHT_Tres_55C_Tres_20C['BHT_Ledaflow (C)']

# Plotting

# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
# #fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# # First subplot: Line plot
# ax1.plot(time_BHP / (24*3600), BHP_Ledaflow / 1e2, color='blue', label='Ledaflow')
# ax1.plot(time_BHP / (24*3600), BHP_GEM / 1e2, color='red', label='GEM', marker='s', markersize=3, linestyle='None')
# #ax1.set_title('Line Plot of Third Column vs Time')
# #ax1.set_xlabel('Time (days)')
# ax1.set_ylabel('BHP (bar)')
# ax1.set_yticks([196, 199, 202, 205, 208, 211])
# ax1.legend()
# ax1.grid(True)

# # Second subplot: Scatter plot
# ax2.plot(time_BHP / (24*3600), BHT_Ledaflow, color='blue', label='Ledaflow')
# ax2.plot(time_BHP / (24*3600), BHT_GEM, color='red', label='GEM', marker='s', markersize=3, linestyle='None')
# ax2.axhline(y=T_SL.iloc[-1], xmin=0, xmax=time_BHT.iloc[-1] / (24*3600), linewidth=2, color='black', linestyle='--')

# #Arrow to indicate
# ax2.annotate(
# # Label and coordinate
# r'$\mathrm{CO}_2\;T_{\mathrm{crit}}\;(\mathrm{NIST}) = 30.9782\,^\circ\mathrm{C}$', xy=(12, 31.2), xytext=(12, 43),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

# #ax2.set_title('Scatter Plot of Third Column vs Time')
# ax2.set_xlabel('Time (days)')
# ax2.set_ylabel(r'BHT ($^\circ$C)')
# ax2.set_yticks([14, 23, 32, 41, 50, 59, 68])
# ax2.legend()
# ax2.grid(True)

# # Adjust layout to prevent overlap and improve appearance
# plt.tight_layout()

# # Save the figure
# plt.savefig("BHP_BHT_results_CO2LINK.png", dpi=300, bbox_inches='tight')

# # Display the plots
# plt.show()

#----------------------------------------------------------------------------

# Plotting

plt.plot(time_WHP / (30 * 24 * 3600), WHP_Ledaflow, color='blue', label='CO2LINK')
plt.plot(twh / (30 * 24 * 3600), pwh / 1e5, color='red', label='T2WELL/ECO2N')
plt.xlabel('Time (months)')
plt.ylabel('WHP (bar)')
plt.yticks([9, 22, 35, 48, 61, 74, 87, 100, 113, 126])
#plt.ylim(0, 100)
plt.grid(True)
plt.legend(loc='best')
plt.tight_layout()
save_figure_png_and_eps(plt.gcf(), "WHP_comparison_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)
plt.show()

# Add horizontal line for CO2 critical density
T_critical_K = CP.PropsSI('Tcrit', 'CO2')  # Critical temperature in K
P_critical_Pa = CP.PropsSI('Pcrit', 'CO2')  # Critical pressure in Pa
rho_critical = CP.PropsSI('D', 'T', T_critical_K, 'P', P_critical_Pa, 'CO2')  # Critical density in kg/m³

plt.plot(time_WH_Density / (30 * 24 * 3600), WH_Density_Ledaflow, color='blue', label='CO2LINK')
plt.plot(twh / (30 * 24 * 3600), dgwh, color='red', label='T2WELL/ECO2N')
plt.axhline(y=rho_critical, color='green', linestyle='--', linewidth=2, label=f'CO₂ Critical Density ({rho_critical:.1f} kg/m³)')
plt.xlabel('Time (months)')
plt.ylabel(r'WH Density (kg/m$^3$)')
#plt.yticks([2, 112, 180, 248, 316, 384, 452, 520, 588, 656, 724, 770])
plt.ylim(0, 800)
plt.grid(True)
plt.legend(loc='best')
plt.tight_layout()
save_figure_png_and_eps(plt.gcf(), "WH_Density_comparison_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)
plt.show()

plt.plot(time_WHP / (30 * 24 * 3600), WH_Density_CoolProp, color='blue', label='CO2LINK')
plt.plot(twh / (30 * 24 * 3600), dgwh, color='red', label='T2WELL/ECO2N')
plt.axhline(y=rho_critical, color='green', linestyle='--', linewidth=2, label=f'CO₂ Critical Density ({rho_critical:.1f} kg/m³)')
plt.xlabel('Time (months)')
plt.ylabel(r'WH Density (kg/m$^3$)')
#plt.yticks([2, 112, 180, 248, 316, 384, 452, 520, 588, 656, 724, 770])
plt.ylim(0, 800)
plt.grid(True)
plt.legend(loc='best')
plt.tight_layout()
save_figure_png_and_eps(plt.gcf(), "WH_Density_comparison_CoolProp_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)
plt.show()

#-----------------------------------

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(3, 4))
#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# First subplot: Line plot
ax1.plot(time_WHP / (24 * 3600), WHP_Ledaflow, color='blue', label='CO2LINK')
ax1.plot(twh / (24 * 3600), pwh * 0.00001, color='red', label='T2WELL/ECO2N')

ax1.set_ylabel('WHP (bar)', fontsize=8)
ax1.set_yticks([9, 22, 35, 48, 61, 74, 87, 100, 113, 126])
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.tick_params(axis='both', which='minor', labelsize=6)
ax1.legend(fontsize=5,loc='best')
ax1.grid(True)

ax2.plot(time_WHP / (24 * 3600), WH_Density_CoolProp, color='blue', label='CO2LINK')
ax2.plot(twh / (24 * 3600), dgwh, color='red', label='T2WELL/ECO2N')
ax2.set_xlabel('Time (days)', fontsize=8)
ax2.set_ylabel(r'WH Density (kg/m$^3$)', fontsize=8)
#ax2.set_yticks([44, 112, 180, 248, 316, 384, 452, 520, 588, 656, 724, 792])
ax2.set_ylim(0, 800)
ax2.tick_params(axis='both', which='major', labelsize=6)
ax2.tick_params(axis='both', which='minor', labelsize=6)
ax2.legend(fontsize=5, loc='best')
ax2.grid(True)

# Adjust layout to prevent overlap and improve appearance
#plt.tight_layout()

# Save the figure
save_figure_png_and_eps(plt.gcf(), "Wellhead_comparison_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)

# Display the plots
plt.show()

#-----------------------------------

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# First subplot: Line plot
#ax1.plot(time_BHP * 30 / 365, BHP_Ledaflow / 1e2, color='blue', label='Ledaflow')
#ax1.plot(time_BHP * 30 / 365, BHP_GEM / 1e2, color='red', label='GEM', marker='s', markersize=3, linestyle='None')
ax1.plot(time_BHP / (24 * 3600), BHP_Ledaflow / 1e2, color='blue', label='CO2LINK')
ax1.plot(tbh_5 / (24 * 3600), pbh_5 * 0.00001, color='red', label='T2WELL/ECO2N')

#Arrow to indicate
ax1.annotate(
# Label and coordinate
r'$\mathrm{Initial}\; P_{res}\; =\; 155\; \mathrm{bar}$', xy=(0.05, 155), xytext=(10, 155),
horizontalalignment="center",
# Custom arrow
arrowprops=dict(arrowstyle='->',lw=1)
)

#ax1.set_title('Line Plot of Third Column vs Time')
#ax1.set_xlabel('Time (days)')
ax1.set_ylabel('BHP (bar)')
ax1.set_yticks([140, 144, 148, 152, 156, 160, 164])
ax1.legend(loc='best')
ax1.grid(True)

# Second subplot: Scatter plot
# ax2.plot(time_BHP[:20] / (24*60), BHT_Ledaflow[:20], color='blue', label='Ledaflow')
# ax2.plot(time_BHP[:20] / (24*60), BHT_GEM[:20], color='red', label='GEM', marker='s', markersize=3, linestyle='None')
ax2.plot(time_BHP / (24 * 3600), BHT_Ledaflow, color='blue', label='CO2LINK')
ax2.plot(tbh_5 / (24 * 3600), tempbh_5, color='red', label='T2WELL/ECO2N')
#ax2.axhline(y=T_SL.iloc[-1], xmin=0, xmax=time_BHT.iloc[-1] / (24*3600), linewidth=2, color='black', linestyle='--')

# #Arrow to indicate
# ax2.annotate(
# # Label and coordinate
# r'$\mathrm{CO}_2\;T_{\mathrm{crit}}\;(\mathrm{NIST}) = 30.9782\,^\circ\mathrm{C}$', xy=(5, 31.2), xytext=(5, 40),
# horizontalalignment="center",
# # Custom arrow
# arrowprops=dict(arrowstyle='->',lw=1)
# )

#Arrow to indicate
ax2.annotate(
# Label and coordinate
r'$\mathrm{Initial}\; T_{res}\; =\; 80\,^\circ\mathrm{C}$', xy=(0.05, 80), xytext=(10, 80),
horizontalalignment="center",
# Custom arrow
arrowprops=dict(arrowstyle='->',lw=1)
)

#ax2.set_title('Scatter Plot of Third Column vs Time')
ax2.set_xlabel('Time (days)')
ax2.set_ylabel(r'BHT ($^\circ$C)')
ax2.set_yticks([45, 50, 55, 60, 65, 70, 75, 80, 85, 90])
ax2.legend(loc='best')
ax2.grid(True)

# Adjust layout to prevent overlap and improve appearance
plt.tight_layout()

# Save the figure
save_figure_png_and_eps(plt.gcf(), "BHP_BHT_comparison_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)

# Display the plots
plt.show()

#-----------------------------------

# Velocity right above the perforations

plt.plot(time_above_perforations_velocity / (30 * 24 * 3600), abs(above_perforations_velocity_Ledaflow), color='blue', label='CO2LINK')
plt.plot(x_above_perforations / (30 * 24 * 3600), abs(y_above_perforations_umix), color='red', label='T2WELL')
#plt.plot(twh / (30 * 24 * 3600), pwh / 1e5, color='red', label='T2WELL/ECO2N')
plt.xlabel('Time (months)')
plt.ylabel('Velocity (m/s)')
#plt.yticks([50, 84, 118, 152, 186, 220, 254, 288, 322, 356, 390])
plt.ylim([-1, 9])
plt.grid(True)
plt.legend(loc='best')
plt.tight_layout()
save_figure_png_and_eps(plt.gcf(), "BH_Velocity_comparison_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)
plt.show()

#----------------------------------------------------------------------------

plt.figure(figsize=(3, 2))
plt.plot(time_BHP / (24 * 3600), BHP_Ledaflow / 1e2, color='blue', label='CO2LINK')
plt.plot(tbh_5 / (24 * 3600), pbh_5 * 0.00001, color='red', label='T2WELL/ECO2N')
#plt.plot(twh / (30 * 24 * 3600), pwh / 1e5, color='red', label='T2WELL/ECO2N')
plt.xlabel('Time (days)', fontsize=7)
plt.ylabel('BHP (bar)', fontsize=7)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.yticks([140, 144, 148, 152, 156, 160, 164])
plt.grid(True)
plt.legend(fontsize=5, loc='best')
plt.tight_layout()
save_figure_png_and_eps(plt.gcf(), "BHP_comparison_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)
plt.show()

#-----------------------------------

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(3, 4))
#fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

# First subplot: Line plot
ax1.plot(time_BHP / (24 * 3600), BHT_Ledaflow, color='blue', label='CO2LINK')
ax1.plot(tbh_5 / (24 * 3600), tempbh_5, color='red', label='T2WELL/ECO2N')

ax1.set_ylabel(r'BHT ($^\circ$C)', fontsize=8)
ax1.set_yticks([45, 50, 55, 60, 65, 70, 75, 80, 85, 90])
ax1.tick_params(axis='both', which='major', labelsize=6)
ax1.tick_params(axis='both', which='minor', labelsize=6)
ax1.legend(fontsize=5, loc='best')
ax1.grid(True)

# Second subplot: Scatter plot
ax2.plot(time_WHP / (24 * 3600), abs(above_perforations_velocity_Ledaflow), color='blue', label='CO2LINK')
ax2.plot(x_above_perforations / (24 * 3600), abs(y_above_perforations_umix), color='red', label='T2WELL/ECO2N')

#Arrow to indicate
ax2.annotate(
# Label and coordinate
r'$\mathrm{Initial}\; T_{res}\; =\; 80\,^\circ\mathrm{C}$', xy=(0.05, 80), xytext=(10, 80),
horizontalalignment="center",
# Custom arrow
arrowprops=dict(arrowstyle='->',lw=1)
)

#ax2.set_title('Scatter Plot of Third Column vs Time')
ax2.set_xlabel('Time (days)', fontsize=8)
ax2.set_ylabel('Velocity (m/s)', fontsize=8)
ax2.tick_params(axis='both', which='major', labelsize=6)
ax2.tick_params(axis='both', which='minor', labelsize=6)
#ax2.set_yticks([40, 45, 50, 55, 60, 65, 70, 75, 80, 85])
ax2.set_ylim([-1, 9])
ax2.legend(fontsize=5, loc='best')
ax2.grid(True)

# Adjust layout to prevent overlap and improve appearance
plt.tight_layout()

# Save the figure
save_figure_png_and_eps(plt.gcf(), "Bottomhole_comparison_T2WELL_CO2LINK_no_caprock_no_salinity_final.png", dpi=300)

# Display the plots
plt.show()