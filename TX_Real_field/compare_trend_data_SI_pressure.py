import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Rectangle, ConnectionPatch
from CoolProp.CoolProp import PropsSI

# Import CoolProp for critical point
try:
    CO2_AVAILABLE = True
    # Get CO2 critical point
    T_CRIT_K = PropsSI('Tcrit', 'CO2')  # Critical temperature in Kelvin
    P_CRIT_PA = PropsSI('Pcrit', 'CO2')  # Critical pressure in Pascal
    T_CRIT_C = T_CRIT_K - 273.15  # Convert to Celsius
    P_CRIT_BAR = P_CRIT_PA / 1e5  # Convert to bar
    print(f"CO2 Critical Point: {T_CRIT_C:.2f}°C, {P_CRIT_BAR:.2f} bar")
except ImportError:
    CO2_AVAILABLE = False
    # Fallback values for CO2 critical point
    T_CRIT_C = 31.0  # °C
    P_CRIT_BAR = 73.8  # bar
    print("CoolProp not available, using fallback critical point values")

# Unit conversions
SEC_TO_DAYS = 1.0 / 86400.0

def save_figure_png_and_eps(fig, filename, dpi=300):
    """Save a Matplotlib figure in both PNG and EPS formats.

    The PNG file keeps the original filename, while the EPS file uses the
    same base name with a .eps extension.
    """
    png_path = Path(filename)
    eps_path = png_path.with_suffix(".eps")

    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(eps_path, format="eps", dpi=dpi, bbox_inches="tight")

    print(f"Saved: {png_path}")
    print(f"Saved: {eps_path}")
    
def plot_with_inset_comparison(liquid_df, sc_df, x_col, y_col_liquid, y_col_sc, 
                              ylabel, filename, ymin_offset=0.0, ymax_offset=0.0,
                              tmax_s_zoom=1000, rect_tmax_s=10000, inset_bbox=(-0.1, -0.15, 1.0, 1.0)):
    """Create comparison plot with zoom-in inset following WH_Results.py pattern."""
    
    tmax_days_rect = rect_tmax_s * SEC_TO_DAYS
    fig, ax = plt.subplots()

    # Plot both datasets on main plot
    ax.plot(liquid_df[x_col], liquid_df[y_col_liquid], 
            label='Liquid CO2', linewidth=2, marker='o', markersize=3, alpha=0.8)
    ax.plot(sc_df[x_col], sc_df[y_col_sc], 
            label='SC CO2', linewidth=2, marker='s', markersize=3, alpha=0.8)
    
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.legend()

    # Rectangle: x=[x0, x0+width], y=[min, max]
    all_values = pd.concat([liquid_df[y_col_liquid], sc_df[y_col_sc]])
    ymin = all_values.min() + ymin_offset
    ymax = all_values.max() + ymax_offset

    rect = Rectangle(
        (-0.05, ymin),
        tmax_days_rect,
        ymax - ymin,
        fill=False,
        linewidth=1.5
    )
    ax.add_patch(rect)

    # Inset
    axins = inset_axes(
        ax,
        width="45%",
        height="45%",
        loc="upper right",
        bbox_to_anchor=inset_bbox,
        bbox_transform=ax.transAxes,
        borderpad=1.2
    )

    # Zoom data (first 1500 seconds)
    liquid_zoom = liquid_df[liquid_df["Time [s]"] <= tmax_s_zoom]
    sc_zoom = sc_df[sc_df["Time [s]"] <= tmax_s_zoom]
    
    axins.plot(liquid_zoom["Time [s]"], liquid_zoom[y_col_liquid], 
              label='Liquid CO2', linewidth=2, marker='o', markersize=3, alpha=0.8)
    axins.plot(sc_zoom["Time [s]"], sc_zoom[y_col_sc], 
              label='SC CO2', linewidth=2, marker='s', markersize=3, alpha=0.8)
    
    axins.set_xlabel("Time [s]", fontsize=8)
    axins.set_ylabel(ylabel, fontsize=8)
    axins.grid(True)
    axins.tick_params(labelsize=8)
    axins.legend(fontsize=8)

    # Connector lines: rectangle top/bottom -> inset top/bottom
    x_left = rect.get_x()
    x_right = x_left + rect.get_width()
    y_bot = rect.get_y()
    y_top = y_bot + rect.get_height()

    x_from = x_right  # connect from right edge of rectangle

    con_top = ConnectionPatch(
        xyA=(0.05, 1.0), coordsA=axins.transAxes,
        xyB=(x_from, y_top), coordsB=ax.transData,
        linewidth=1.2, color="0.35"
    )
    con_bot = ConnectionPatch(
        xyA=(0.05, 0.0), coordsA=axins.transAxes,
        xyB=(x_from, y_bot), coordsB=ax.transData,
        linewidth=1.2, color="0.35"
    )

    ax.add_artist(con_top)
    ax.add_artist(con_bot)

    plt.tight_layout()
    save_figure_png_and_eps(fig, filename, dpi=300)
    plt.show()

def plot_single_dataset(df, x_col, y_col, ylabel, filename):
    """Create plot for a single dataset when comparison is not possible."""
    fig, ax = plt.subplots()
    
    ax.plot(df[x_col], df[y_col], 
            label='Liquid CO2', linewidth=2, marker='o', markersize=3, alpha=0.8)
    
    ax.set_xlabel("Time [days]")
    ax.set_ylabel(ylabel)
    ax.grid(True)
    ax.legend()
    
    plt.tight_layout()
    save_figure_png_and_eps(fig, filename, dpi=300)
    plt.show()
    print(f"Saved: {filename}") 

def load_wh_files(folder_path):
    """Load and process wellhead files following WH_Results.py pattern."""
    folder = Path(folder_path)
    
    # Load WH files with specific column selection and skip rows
    try:
        df_mass = pd.read_csv(folder / "WH_MassFlowRate.csv", usecols=[0, 1, 2, 3, 4, 5, 6, 7], skiprows=2, header=None)
        df_whp  = pd.read_csv(folder / "WHP.csv", usecols=[0, 1], skiprows=2, header=None)
        df_whd  = pd.read_csv(folder / "WH_Density.csv", usecols=[0, 1, 2, 3, 4, 5], skiprows=2, header=None)
        
        # Load WH_Viscosity if available
        df_whv = None
        if (folder / "WH_Viscosity.csv").exists():
            df_whv = pd.read_csv(folder / "WH_Viscosity.csv", usecols=[0, 1, 2, 3, 4], skiprows=2, header=None)
        
        # Set column names
        df_mass.columns = ["Time [s]", "Mass Rate - Bubbles [kg/s]", "Mass Rate - Continuous Gas [kg/s]", "Mass Rate - Continuous Liquid [kg/s]", "Mass Rate - Droplets [kg/s]", "Mass Rate - Total [kg/s]", "Mass Rate - Total Gas [kg/s]", "Mass Rate - Total Liquid [kg/s]"]
        df_whp.columns  = ["Time [s]", "Pressure [bar]"]
        df_whd.columns  = ["Time [s]", "Density - Bubbles [kg/m3]", "Density - Continuous Gas [kg/m3]", "Density - Continuous Liquid [kg/m3]", "Density - Droplets [kg/m3]", "Density - Mixture [kg/m3]"]
        
        # Set WH_Viscosity column names if loaded
        if df_whv is not None:
            df_whv.columns  = ["Time [s]", "Viscosity - Bubbles [Pa.s]", "Viscosity - Continuous Gas [Pa.s]", "Viscosity - Continuous Liquid [Pa.s]", "Viscosity - Droplets [Pa.s]"]
        
        # Switch sign of mass rate columns
        for col in df_mass.columns:
            if "Mass Rate" in col:
                df_mass[col] = -df_mass[col]
        
        # Apply conversions
        df_whp["Time [days]"] = df_whp["Time [s]"] * SEC_TO_DAYS
        df_mass["Time [days]"] = df_mass["Time [s]"] * SEC_TO_DAYS
        df_whd["Time [days]"] = df_whd["Time [s]"] * SEC_TO_DAYS
        if df_whv is not None:
            df_whv["Time [days]"] = df_whv["Time [s]"] * SEC_TO_DAYS
        
        # Load WHT if available
        df_wht = None
        if (folder / "WHT.csv").exists():
            df_wht = pd.read_csv(folder / "WHT.csv", usecols=[0, 1, 3, 4], skiprows=2, header=None)
            df_wht.columns = ["Time [s]", "Temperature - Average[C]", "Temperature - Gas[C]", "Temperature - Liquid[C]"]
            df_wht["Time [days]"] = df_wht["Time [s]"] * SEC_TO_DAYS
            # Add Fahrenheit column for WHT
            df_wht["Temperature [°F]"] = df_wht["Temperature - Average[C]"]
        
        return_dict = {
            'WHP': df_whp,
            'WH_MassFlowRate': df_mass,
            'WH_Density': df_whd,
            'WHT': df_wht
        }
        
        # Add WH_Viscosity if loaded
        if df_whv is not None:
            return_dict['WH_Viscosity'] = df_whv
            
        return return_dict
    
    except FileNotFoundError as e:
        print(f"Missing WH file in {folder_path}: {e}")
        return {}

def load_bh_files(folder_path):
    """Load and process bottomhole files following BH_Results.py pattern."""
    folder = Path(folder_path)
    
    try:
        df_bhp = pd.read_csv(folder / "BHP.csv", usecols=[0, 1], skiprows=2, header=None)
        df_bhp.columns = ["Time [s]", "Pressure [kPa]"]
        df_bhp["Pressure [bar]"] = df_bhp["Pressure [kPa]"] / 100.0  # Convert kPa to bar
        df_bhp["Time [days]"] = df_bhp["Time [s]"] * SEC_TO_DAYS
        
        result = {'BHP': df_bhp}
        
        # Load BHT if available
        if (folder / "BHT.csv").exists():
            df_bht = pd.read_csv(folder / "BHT.csv", usecols=[0, 1], skiprows=2, header=None)
            df_bht.columns = ["Time [s]", "Temperature [C]"]
            df_bht["Time [days]"] = df_bht["Time [s]"] * SEC_TO_DAYS
            result['BHT'] = df_bht
            
            # Calculate BH Density and BH Viscosity using CoolProp if both BHP and BHT are available
            if CO2_AVAILABLE:
                print(f"BHT data loaded with {len(df_bht)} rows")
                print(f"BHP data loaded with {len(df_bhp)} rows")
                try:
                    # Merge BHP and BHT data on time (remove duplicate Time [days] column first)
                    df_bht_merge = df_bht.drop('Time [days]', axis=1)  # Remove duplicate column
                    df_merged = pd.merge(df_bhp, df_bht_merge, on='Time [s]', how='inner')
                    print(f"Merged data has {len(df_merged)} rows")
                    
                    # Calculate density and viscosity using CoolProp
                    bh_density = []
                    bh_viscosity = []
                    
                    for _, row in df_merged.iterrows():
                        temp_c = row['Temperature [C]']
                        pressure_kpa = row['Pressure [kPa]']  # Use original kPa value
                        
                        # CoolProp calculations (temperature in K, pressure in Pa)
                        density = PropsSI('D', 'T', temp_c + 273.15, 'P', pressure_kpa * 1000, 'CO2')
                        viscosity = PropsSI('V', 'T', temp_c + 273.15, 'P', pressure_kpa * 1000, 'CO2')
                        
                        bh_density.append(density)
                        bh_viscosity.append(viscosity)
                    
                    # Create dataframes for BH properties
                    df_bh_density = pd.DataFrame({
                        'Time [s]': df_merged['Time [s]'],
                        'Time [days]': df_merged['Time [days]'],
                        'Density [kg/m3]': bh_density
                    })
                    
                    df_bh_viscosity = pd.DataFrame({
                        'Time [s]': df_merged['Time [s]'],
                        'Time [days]': df_merged['Time [days]'],
                        'Viscosity [Pa.s]': bh_viscosity
                    })
                    
                    result['BH_Density'] = df_bh_density
                    result['BH_Viscosity'] = df_bh_viscosity
                    
                    print(f"Calculated BH Density and BH Viscosity using CoolProp for {len(bh_density)} data points")
                    print(f"BH_Density range: {min(bh_density):.1f} to {max(bh_density):.1f} kg/m3")
                    print(f"BH_Viscosity range: {min(bh_viscosity):.2e} to {max(bh_viscosity):.2e} Pa.s")
                    
                except Exception as e:
                    print(f"Error calculating BH properties with CoolProp: {e}")
            else:
                print("CoolProp not available, skipping BH Density and BH Viscosity calculations")
            
        return result
        
    except FileNotFoundError as e:
        print(f"Missing BH file in {folder_path}: {e}")
        return {}

def determine_viscosity_column(df_temp, df_pressure, df_viscosity):
    """Determine which viscosity column to use based on CO2 critical point for each data point."""
    # Ensure dataframes are aligned by index (time)
    df_combined = pd.DataFrame({
        'Temperature - Average[C]': df_temp['Temperature - Average[C]'],
        'Pressure [bar]': df_pressure['Pressure [bar]'],
        'Viscosity - Continuous Liquid [Pa.s]': df_viscosity['Viscosity - Continuous Liquid [Pa.s]'],
        'Viscosity - Continuous Gas [Pa.s]': df_viscosity['Viscosity - Continuous Gas [Pa.s]']
    }).dropna()

    selected_viscosity = []
    for index, row in df_combined.iterrows():
        temp_c = row['Temperature - Average[C]']
        pressure_bar = row['Pressure [bar]']

        # Since pressure is always above critical, use temperature to decide
        if temp_c < T_CRIT_C:
            selected_viscosity.append(row['Viscosity - Continuous Liquid [Pa.s]'])
        else:
            selected_viscosity.append(row['Viscosity - Continuous Gas [Pa.s]'])
            
    print(f"Critical point: Tcrit = {T_CRIT_C:.2f}°C, Pcrit = {P_CRIT_BAR:.2f} bar")
    print(f"Dynamically selected viscosity column based on temperature relative to critical point.")
    return pd.Series(selected_viscosity, index=df_combined.index)

def determine_density_column(df_temp, df_pressure, df_density):
    """Determine which density column to use based on CO2 critical point for each data point."""
    # Ensure dataframes are aligned by index (time)
    df_combined = pd.DataFrame({
        'Temperature - Average[C]': df_temp['Temperature - Average[C]'],
        'Pressure [bar]': df_pressure['Pressure [bar]'],
        'Density - Continuous Liquid [kg/m3]': df_density['Density - Continuous Liquid [kg/m3]'],
        'Density - Continuous Gas [kg/m3]': df_density['Density - Continuous Gas [kg/m3]']
    }).dropna()

    selected_density = []
    transition_count = 0
    for index, row in df_combined.iterrows():
        temp_c = row['Temperature - Average[C]']
        pressure_bar = row['Pressure [bar]']
        liquid_density = row['Density - Continuous Liquid [kg/m3]']
        gas_density = row['Density - Continuous Gas [kg/m3]']

        # Since pressure is always above critical, use temperature to decide
        if temp_c < T_CRIT_C:
            selected_density.append(liquid_density)
            if transition_count < 5:  # Print first few transitions
                print(f"Time {index}: T={temp_c:.1f}°C < Tcrit, using LIQUID density={liquid_density:.1f}")
                transition_count += 1
        else:
            selected_density.append(gas_density)
            if transition_count == 0:  # Print first point above critical
                print(f"Time {index}: T={temp_c:.1f}°C >= Tcrit, using GAS density={gas_density:.1f}")
            
    print(f"Critical point: Tcrit = {T_CRIT_C:.2f}°C, Pcrit = {P_CRIT_BAR:.2f} bar")
    print(f"Temperature range: {df_combined['Temperature - Average[C]'].min():.1f} to {df_combined['Temperature - Average[C]'].max():.1f}°C")
    print(f"Liquid density range: {df_combined['Density - Continuous Liquid [kg/m3]'].min():.1f} to {df_combined['Density - Continuous Liquid [kg/m3]'].max():.1f} kg/m3")
    print(f"Gas density range: {df_combined['Density - Continuous Gas [kg/m3]'].min():.1f} to {df_combined['Density - Continuous Gas [kg/m3]'].max():.1f} kg/m3")
    return pd.Series(selected_density, index=df_combined.index)

def plot_comparison(liquid_data, sc_data, parameter_name):
    """Create comparison plot with zoom-in inset for a specific parameter."""
    
    # Check if parameter exists in at least one dataset
    liquid_exists = parameter_name in liquid_data and liquid_data[parameter_name] is not None
    sc_exists = parameter_name in sc_data and sc_data[parameter_name] is not None
    
    if liquid_exists or sc_exists:
        liquid_df = liquid_data[parameter_name] if liquid_exists else None
        sc_df = sc_data[parameter_name] if sc_exists else None
        
        # Use Time [days] for x-axis
        x_col = "Time [days]"
        
        # Determine which column to plot based on parameter
        if parameter_name == 'WHP':
            y_col_liquid = 'Pressure [bar]'
            y_col_sc = 'Pressure [bar]'
            ylabel = 'WHP [bar]'
            ymin_offset = -4
            ymax_offset = -25
            inset_bbox = (-0.27, -0.13, 1.0, 1.0)
        elif parameter_name == 'BHP':
            y_col_liquid = 'Pressure [bar]'
            y_col_sc = 'Pressure [bar]'
            ylabel = 'BHP [bar]'
            ymin_offset = -4
            ymax_offset = -15
            inset_bbox = (-0.29, -0.25, 1.0, 1.0)
        elif parameter_name == 'WH_MassFlowRate':
            # Use Total [kg/s] for both liquid and SC cases
            y_col_liquid = 'Mass Rate - Total [kg/s]'
            y_col_sc = 'Mass Rate - Total [kg/s]'
            ylabel = 'WH Mass Flow Rate [kg/s]'
            ymin_offset = -1
            ymax_offset = -5
            inset_bbox = (-0.25, -0.23, 1.0, 1.0)
        elif parameter_name == 'WH_Density':
            # Create combined density series based on critical point
            if 'WHT' in liquid_data and liquid_data['WHT'] is not None and 'WHP' in liquid_data and liquid_data['WH_Density'] is not None:
                density_series_liquid = determine_density_column(
                    liquid_data['WHT'], liquid_data['WHP'], liquid_data['WH_Density']
                )
                # Create a temporary dataframe with the combined series
                liquid_df = pd.DataFrame({
                    'Time [days]': liquid_data['WHT']['Time [days]'],
                    'Time [s]': liquid_data['WHT']['Time [s]'],
                    'Density [kg/m3]': density_series_liquid
                }).dropna()
                y_col_liquid = 'Density [kg/m3]'
            else:
                y_col_liquid = "Density - Continuous Liquid [kg/m3]"
            
            if 'WHT' in sc_data and sc_data['WHT'] is not None and 'WHP' in sc_data and sc_data['WH_Density'] is not None:
                density_series_sc = determine_density_column(
                    sc_data['WHT'], sc_data['WHP'], sc_data['WH_Density']
                )
                # Create a temporary dataframe with the combined series
                sc_df = pd.DataFrame({
                    'Time [days]': sc_data['WHT']['Time [days]'],
                    'Time [s]': sc_data['WHT']['Time [s]'],
                    'Density [kg/m3]': density_series_sc
                }).dropna()
                y_col_sc = 'Density [kg/m3]'
            else:
                y_col_sc = "Density - Continuous Gas [kg/m3]"
            
            ylabel = 'WH Density [kg/m³]'
            ymin_offset = -30
            ymax_offset = 10
            inset_bbox = (-0.30, -0.27, 1.0, 1.0)
        elif parameter_name == 'WH_Viscosity':
            # Create combined viscosity series based on critical point
            if 'WHT' in liquid_data and liquid_data['WHT'] is not None and 'WHP' in liquid_data and liquid_data['WH_Viscosity'] is not None:
                viscosity_series_liquid = determine_viscosity_column(
                    liquid_data['WHT'], liquid_data['WHP'], liquid_data['WH_Viscosity']
                )
                # Create a temporary dataframe with the combined series
                liquid_df = pd.DataFrame({
                    'Time [days]': liquid_data['WHT']['Time [days]'],
                    'Time [s]': liquid_data['WHT']['Time [s]'],
                    'Viscosity [Pa.s]': viscosity_series_liquid
                }).dropna()
                y_col_liquid = 'Viscosity [Pa.s]'
            else:
                y_col_liquid = "Viscosity - Continuous Liquid [Pa.s]"
            
            if 'WHT' in sc_data and sc_data['WHT'] is not None and 'WHP' in sc_data and sc_data['WH_Viscosity'] is not None:
                viscosity_series_sc = determine_viscosity_column(
                    sc_data['WHT'], sc_data['WHP'], sc_data['WH_Viscosity']
                )
                # Create a temporary dataframe with the combined series
                sc_df = pd.DataFrame({
                    'Time [days]': sc_data['WHT']['Time [days]'],
                    'Time [s]': sc_data['WHT']['Time [s]'],
                    'Viscosity [Pa.s]': viscosity_series_sc
                }).dropna()
                y_col_sc = 'Viscosity [Pa.s]'
            else:
                y_col_sc = "Viscosity - Continuous Gas [Pa.s]"
            
            ylabel = 'WH Viscosity [Pa.s]'
            ymin_offset = -0.000003
            ymax_offset = -0.000001
            inset_bbox = (-0.30, -0.19, 1.0, 1.0)
        elif parameter_name == 'WHT':
            y_col_liquid = 'Temperature [°F]'
            y_col_sc = 'Temperature [°F]'
            ylabel = 'WHT [°F]'
            ymin_offset = -1.1
            ymax_offset = 0.9
            inset_bbox = (-0.29, -0.11, 1.0, 1.0)
        elif parameter_name == 'BHT':
            y_col_liquid = 'Temperature [C]'
            y_col_sc = 'Temperature [C]'
            ylabel = 'BHT [°C]'
            ymin_offset = 7.8
            ymax_offset = 0.3
            inset_bbox = (-0.27, -0.21, 1.0, 1.0)
        elif parameter_name == 'BH_Density':
            y_col_liquid = 'Density [kg/m3]'
            y_col_sc = 'Density [kg/m3]'
            ylabel = 'BH Density [kg/m³]'
            ymin_offset = -6
            ymax_offset = -56
            inset_bbox = (-0.17, -0.13, 1.0, 1.0)
        elif parameter_name == 'BH_Viscosity':
            y_col_liquid = 'Viscosity [Pa.s]'
            y_col_sc = 'Viscosity [Pa.s]'
            ylabel = 'BH Viscosity [Pa.s]'
            ymin_offset = -0.000001
            #ymax_offset = 0.0000005
            ymax_offset = -0.000012
            inset_bbox = (-0.19, -0.16, 1.0, 1.0)
        else:
            # Fallback: find first non-time column
            liquid_data_cols = [col for col in liquid_df.columns if "Time" not in col]
            sc_data_cols = [col for col in sc_df.columns if "Time" not in col]
            y_col_liquid = liquid_data_cols[0] if liquid_data_cols else None
            y_col_sc = sc_data_cols[0] if sc_data_cols else None
            ylabel = parameter_name
            ymin_offset = 0.0
            ymax_offset = 0.0
            inset_bbox = (-0.1, -0.15, 1.0, 1.0)
        
        # Check if columns exist before plotting
        if liquid_df is not None and y_col_liquid and y_col_liquid in liquid_df.columns:
            filename = f'Comparison_{parameter_name}.png'
            
            # If SC data exists and has column, plot comparison
            if sc_df is not None and y_col_sc and y_col_sc in sc_df.columns:
                plot_with_inset_comparison(
                    liquid_df, sc_df, x_col, y_col_liquid, y_col_sc, 
                    ylabel, filename, ymin_offset, ymax_offset, 
                    inset_bbox=inset_bbox
                )
            else:
                # Plot only available data (single dataset)
                plot_single_dataset(liquid_df, x_col, y_col_liquid, ylabel, filename)
        elif sc_df is not None and y_col_sc and y_col_sc in sc_df.columns:
            # Plot only SC data (single dataset)
            filename = f'Comparison_{parameter_name}.png'
            plot_single_dataset(sc_df, x_col, y_col_sc, ylabel, filename)
        else:
            print(f"Required columns not found for {parameter_name}")
    else:
        print(f"Parameter {parameter_name} not found in any dataset")

def main():
    # Define paths
    liquid_trend_path = r"D:\CO2LINK\Paper\Real_field_model\TX_BrownPelican_field\Pressure_Control\For_Feb25\Liquid_CO2_injection\Results_for_Feb25\Trend"
    sc_trend_path = r"D:\CO2LINK\Paper\Real_field_model\TX_BrownPelican_field\Pressure_Control\For_Feb25\SC_CO2_injection\Results_for_Feb25\Trend"
    
    # Load WH and BH files separately
    print("Loading Liquid CO2 data...")
    liquid_wh_data = load_wh_files(liquid_trend_path)
    liquid_bh_data = load_bh_files(liquid_trend_path)
    print(f"Liquid WH data keys: {list(liquid_wh_data.keys())}")
    print(f"Liquid BH data keys: {list(liquid_bh_data.keys())}")
    
    print("\nLoading SC CO2 data...")
    sc_wh_data = load_wh_files(sc_trend_path)
    sc_bh_data = load_bh_files(sc_trend_path)
    print(f"SC WH data keys: {list(sc_wh_data.keys())}")
    print(f"SC BH data keys: {list(sc_bh_data.keys())}")
    
    # Combine all data
    liquid_data = {**liquid_wh_data, **liquid_bh_data}
    sc_data = {**sc_wh_data, **sc_bh_data}
    
    # Get all parameters (exclude None values)
    liquid_params = {key for key, value in liquid_data.items() if value is not None}
    sc_params = {key for key, value in sc_data.items() if value is not None}
    all_params = liquid_params.union(sc_params)
    
    print(f"\nAll parameters found: {sorted(all_params)}")
    print(f"Liquid CO2 only: {sorted(liquid_params - sc_params)}")
    print(f"SC CO2 only: {sorted(sc_params - liquid_params)}")
    
    # Create comparison plots for all parameters
    for param in sorted(all_params):
        print(f"\nCreating comparison plot for {param}...")
        plot_comparison(liquid_data, sc_data, param)
    
    print("\nComparison complete! Plots saved as PNG files.")

if __name__ == "__main__":
    main()
