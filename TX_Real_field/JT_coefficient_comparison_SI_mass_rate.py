# -*- coding: utf-8 -*-
"""
JT Coefficient Comparison - Valve Outlet Temperature Calculator
Calculates outlet temperature for isenthalpic (Joule-Thomson) expansion of CO2
for both Supercritical and Liquid injection cases.

Created on Tue Feb 24 13:57:58 2026
@author: jpauya1
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from CoolProp.CoolProp import PropsSI, PhaseSI

# --------------------------------------------------
# UNIT CONVERSIONS
# --------------------------------------------------

SEC_TO_DAYS = 1.0 / 86400.0
C_TO_F_SCALE = 9/5
C_TO_F_OFFSET = 32
LBFT3_TO_KGM3 = 1/16.018 
PSI_TO_PA = 6894.76
BAR_TO_PA = 1e5

FLUID = "CO2"


def save_figure_png_and_eps(fig, filename, dpi=300):
    """Save a Matplotlib figure in both PNG and EPS formats."""
    png_path = Path(filename)
    eps_path = png_path.with_suffix(".eps")

    fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(eps_path, format="eps", dpi=dpi, bbox_inches="tight")

    print(f"Saved: {png_path}")
    print(f"Saved: {eps_path}")
    
    
def c_to_k(T_c):
    """Convert Celsius to Kelvin."""
    return T_c + 273.15


def k_to_c(T_k):
    """Convert Kelvin to Celsius."""
    return T_k - 273.15


def calculate_valve_outlet_temperature(P1_data, P2_data, Tin_C, case_name="Case"):
    """
    Calculate valve outlet temperature using isenthalpic expansion.
    
    Parameters:
    -----------
    P1_data : array-like
        Inlet pressure values (in bar)
    P2_data : array-like
        Outlet pressure values (in bar)
    Tin_C : float
        Inlet temperature in Celsius
    case_name : str
        Name of the case (for reporting)
    
    Returns:
    --------
    dict : Contains T2_C (outlet temp in C), phases, qualities, and status
    """
    
    print(f"\n=== Processing {case_name} ===")
    print(f"Inlet temperature: {Tin_C}°C")
    print(f"Number of data points: {len(P1_data)}")
    
    # Convert inlet temperature to Kelvin
    Tin_K = c_to_k(Tin_C)
    
    # Convert pressures to Pascal for CoolProp
    P1_pa = np.array(P1_data) * BAR_TO_PA
    P2_pa = np.array(P2_data) * BAR_TO_PA
    
    # Calculate inlet enthalpy for each point
    print("Calculating inlet enthalpies...")
    h1 = np.zeros(len(P1_pa))
    for i in range(len(P1_pa)):
        try:
            h1[i] = PropsSI("H", "P", P1_pa[i], "T", Tin_K, FLUID)
        except Exception as e:
            print(f"  Warning: Failed to calculate h1 at point {i}: P1={P1_data[i]:.2f} bar, T={Tin_C}°C")
            h1[i] = np.nan
    
    # Calculate outlet state using (P, H) - isenthalpic condition
    print("Calculating outlet states...")
    T2_list = []
    Q_list = []
    phase_list = []
    
    # Ensure we only process the minimum length to avoid index errors
    n_points = min(len(P2_pa), len(h1))
    print(f"Processing {n_points} points (P2: {len(P2_pa)}, h1: {len(h1)})")
    
    for i in range(n_points):
        if i % 1000 == 0:
            print(f"  Processing point {i}/{n_points}")
        
        if np.isnan(h1[i]):
            T2_list.append(np.nan)
            Q_list.append(np.nan)
            phase_list.append("FAILED")
            continue
        
        try:
            # Temperature from isenthalpic condition
            T2 = PropsSI("T", "P", P2_pa[i], "H", h1[i], FLUID)
            T2_list.append(T2)
            
            # Detect phase
            phase = PhaseSI("P", P2_pa[i], "H", h1[i], FLUID)
            phase_list.append(phase)
            
            # Try to get vapor quality (only valid in two-phase)
            try:
                Q = PropsSI("Q", "P", P2_pa[i], "H", h1[i], FLUID)
            except:
                Q = np.nan
            Q_list.append(Q)
            
        except Exception as e:
            print(f"  Warning: Failed at point {i}: P2={P2_data[i]:.2f} bar, h1={h1[i]:.2f} J/kg")
            T2_list.append(np.nan)
            Q_list.append(np.nan)
            phase_list.append("FAILED")
    
    # Convert to arrays
    T2_array = np.array(T2_list)
    Q_array = np.array(Q_list)
    
    # Convert outlet temperature to Celsius
    T2_C = k_to_c(T2_array)
    
    # Print summary
    print(f"\n{case_name} Summary:")
    print(f"  Outlet temperature range: {np.nanmin(T2_C):.2f}°C to {np.nanmax(T2_C):.2f}°C")
    unique_phases = set([p for p in phase_list if p != "FAILED"])
    print(f"  Phases encountered: {unique_phases}")
    
    return {
        'T2_C': T2_C,
        'T2_K': T2_array,
        'Q': Q_array,
        'phases': phase_list,
        'h1': h1,
        'status': 'SUCCESS' if not np.all(np.isnan(T2_C)) else 'FAILED'
    }


def load_pressure_data(file_path, skiprows=2, usecols=None):
    """
    Load pressure data from CSV file.
    
    Parameters:
    -----------
    file_path : str or Path
        Path to the CSV file
    skiprows : int
        Number of rows to skip (header rows)
    usecols : list or None
        Columns to use (default: first two columns)
    
    Returns:
    --------
    tuple : (time_array, pressure_array) where pressure is in bar
    """
    try:
        df = pd.read_csv(file_path, skiprows=skiprows, header=None)
        
        if usecols is None:
            usecols = [0, 1]  # Default: first two columns
        
        # Extract time and pressure
        time = df.iloc[:, usecols[0]].values
        pressure = df.iloc[:, usecols[1]].values
        
        print(f"Loaded {file_path}: {len(pressure)} data points")
        return time, pressure
        
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None


def load_temperature_data(file_path, skiprows=2, usecols=None):
    """
    Load temperature data from CSV file.
    
    Parameters:
    -----------
    file_path : str or Path
        Path to the CSV file
    skiprows : int
        Number of rows to skip (header rows)
    usecols : list or None
        Columns to use (default: first two columns)
    
    Returns:
    --------
    tuple : (time_array, temperature_array) where temperature is in original units
    """
    try:
        df = pd.read_csv(file_path, skiprows=skiprows, header=None)
        
        if usecols is None:
            usecols = [0, 1]  # Default: first two columns
        
        # Extract time and temperature
        time = df.iloc[:, usecols[0]].values
        temperature = df.iloc[:, usecols[1]].values
        
        print(f"Loaded {file_path}: {len(temperature)} data points")
        return time, temperature
        
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None, None



def plot_temperature_comparison(time_sc, T2_sc_C, time_liquid, T2_liquid_C, 
                               time_sc_wht, T_sc_wht_C, time_liquid_wht, T_liquid_wht_C,
                               case_sc_name="SC Case", case_liquid_name="Liquid Case",
                               filename="JT_outlet_temperature_comparison.png"):
    """
    Plot comparison of outlet temperatures for SC and Liquid cases with actual WHT data.
    """
    
    fig = plt.figure()
    
    # Plot JT calculated temperatures with crosses
    plt.plot(time_sc/(24*3600), T2_sc_C, label=f'{case_sc_name} (JT calculated)', 
             color='red')
    plt.plot(time_liquid/(24*3600), T2_liquid_C, label=f'{case_liquid_name} (JT calculated)', 
             color='blue')
    
    # Plot actual WHT data with solid lines
    if time_sc_wht is not None and T_sc_wht_C is not None:
        plt.plot(time_sc_wht/(24*3600), T_sc_wht_C, label=f'{case_sc_name} (WHT)', 
                 color='red', linestyle='dashed')
    
    if time_liquid_wht is not None and T_liquid_wht_C is not None:
        plt.plot(time_liquid_wht/(24*3600), T_liquid_wht_C, label=f'{case_liquid_name} (WHT)',
                 color='blue', linestyle='dashed')
    
    plt.xlabel("Time [days]")
    plt.ylabel("Temperature [°C]")
    plt.title("Valve Outlet Temperature")
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    plt.tight_layout()
    save_figure_png_and_eps(fig, filename, dpi=300)
    plt.show()


def main():
    """
    Main function to calculate and compare valve outlet temperatures
    for SC and Liquid CO2 injection cases.
    """
    
    print("="*70)
    print("JT COEFFICIENT COMPARISON - VALVE OUTLET TEMPERATURE")
    print("Isenthalpic (Joule-Thomson) Expansion Analysis")
    print("="*70)
    
    # --------------------------------------------------
    # DEFINE PATHS
    # --------------------------------------------------
    # Using MassRate_control paths (as in original file)
    liquid_trend_path = r"d:\CO2LINK\Paper\Real_field_model\TX_BrownPelican_field\MassRate_Control\For_Feb25\Liquid_CO2_injection\Results_for_Feb25\Trend"
    sc_trend_path = r"d:\CO2LINK\Paper\Real_field_model\TX_BrownPelican_field\MassRate_Control\For_Feb25\SC_CO2_injection\Results_for_Feb25\Trend"
    
    # --------------------------------------------------
    # USER INPUTS - Inlet Temperatures
    # --------------------------------------------------
    # SC case: 95°F inlet temperature
    # Liquid case: 70°F inlet temperature
    Tin_SC_F = 95.0
    Tin_Liquid_F = 70.0
    
    # Convert to Celsius using existing conversion constants
    Tin_SC_C = (Tin_SC_F - C_TO_F_OFFSET) / C_TO_F_SCALE
    Tin_Liquid_C = (Tin_Liquid_F - C_TO_F_OFFSET) / C_TO_F_SCALE
    
    print(f"\nInlet Conditions:")
    print(f"  SC CO2: {Tin_SC_F}°F ({Tin_SC_C:.1f}°C)")
    print(f"  Liquid CO2: {Tin_Liquid_F}°F ({Tin_Liquid_C:.1f}°C)")
    
    # --------------------------------------------------
    # LOAD PRESSURE DATA
    # --------------------------------------------------
    print("\n" + "-"*50)
    print("LOADING PRESSURE DATA")
    print("-"*50)
    
    # SC Case - Inlet (BHP) and Outlet (WHP) pressures
    print("\nSC Case:")
    P1_file_sc = os.path.join(sc_trend_path, "ValveInletPressure.csv")  # Inlet: Bottom Hole Pressure
    P2_file_sc = os.path.join(sc_trend_path, "ValveOutletPressure.csv")  # Outlet: Well Head Pressure
    
    time_sc_in, P1_sc = load_pressure_data(P1_file_sc)
    time_sc_out, P2_sc = load_pressure_data(P2_file_sc)
    
    # Liquid Case - Inlet (BHP) and Outlet (WHP) pressures
    print("\nLiquid Case:")
    P1_file_liquid = os.path.join(liquid_trend_path, "ValveInletPressure.csv")  # Inlet
    P2_file_liquid = os.path.join(liquid_trend_path, "ValveOutletPressure.csv")  # Outlet
    
    time_liquid_in, P1_liquid = load_pressure_data(P1_file_liquid)
    time_liquid_out, P2_liquid = load_pressure_data(P2_file_liquid)
    
    # Check if data loaded successfully
    if P1_sc is None or P2_sc is None:
        print("\nERROR: Failed to load SC pressure data!")
        return
    
    if P1_liquid is None or P2_liquid is None:
        print("\nERROR: Failed to load Liquid pressure data!")
        return
    
    # --------------------------------------------------
    # LOAD TEMPERATURE DATA (WHT)
    # --------------------------------------------------
    print("\n" + "-"*50)
    print("LOADING TEMPERATURE DATA")
    print("-"*50)
    
    # Load WHT data for SC case
    print("\nSC Case WHT:")
    WHT_file_sc = os.path.join(sc_trend_path, "WHT.csv")
    time_sc_wht, T_sc_wht = load_temperature_data(WHT_file_sc)
    
    # Load WHT data for Liquid case
    print("\nLiquid Case WHT:")
    WHT_file_liquid = os.path.join(liquid_trend_path, "WHT.csv")
    time_liquid_wht, T_liquid_wht = load_temperature_data(WHT_file_liquid)
    
    # WHT.csv exports Celsius, as declared in its header; do not infer units
    # from magnitude (valid Celsius values can exceed 50).
    T_sc_wht_C = T_sc_wht
    T_liquid_wht_C = T_liquid_wht

    # --------------------------------------------------
    # CALCULATE OUTLET TEMPERATURES
    # --------------------------------------------------
    print("\n" + "-"*50)
    print("CALCULATING OUTLET TEMPERATURES")
    print("-"*50)
    
    # SC Case
    results_sc = calculate_valve_outlet_temperature(
        P1_sc, P2_sc, Tin_SC_C, case_name="SC Case"
    )
    
    # Liquid Case
    results_liquid = calculate_valve_outlet_temperature(
        P1_liquid, P2_liquid, Tin_Liquid_C, case_name="Liquid Case"
    )
    
    # --------------------------------------------------
    # SAVE RESULTS
    # --------------------------------------------------
    print("\n" + "-"*50)
    print("SAVING RESULTS")
    print("-"*50)
    
    # Create DataFrames with results
    df_sc = pd.DataFrame({
        'Time [s]': time_sc_out if time_sc_out is not None else range(len(results_sc['T2_C'])),
        'Inlet_Pressure [bar]': P1_sc[:len(results_sc['T2_C'])],
        'Outlet_Pressure [bar]': P2_sc,
        'Inlet_Temperature [C]': Tin_SC_C,
        'Outlet_Temperature [C]': results_sc['T2_C'],
        'Wellhead_Temperature [C]': T_sc_wht_C,
        'Phase': results_sc['phases'],
        'Vapor_Quality': results_sc['Q']
    })
    
    df_liquid = pd.DataFrame({
        'Time [s]': time_liquid_out if time_liquid_out is not None else range(len(results_liquid['T2_C'])),
        'Inlet_Pressure [bar]': P1_liquid[:len(results_liquid['T2_C'])],
        'Outlet_Pressure [bar]': P2_liquid,
        'Inlet_Temperature [C]': Tin_Liquid_C,
        'Outlet_Temperature [C]': results_liquid['T2_C'],
        'Wellhead_Temperature [C]': T_liquid_wht_C,
        'Phase': results_liquid['phases'],
        'Vapor_Quality': results_liquid['Q']
    })
    
    # Save to CSV
    output_sc = "JT_outlet_temperature_SC.csv"
    output_liquid = "JT_outlet_temperature_Liquid.csv"
    
    df_sc.to_csv(output_sc, index=False)
    df_liquid.to_csv(output_liquid, index=False)
    
    print(f"SC results saved to: {output_sc}")
    print(f"Liquid results saved to: {output_liquid}")
    
    # --------------------------------------------------
    # PLOT RESULTS
    # --------------------------------------------------
    print("\n" + "-"*50)
    print("GENERATING PLOTS")
    print("-"*50)
    
    plot_temperature_comparison(
        time_sc_out if time_sc_out is not None else range(len(results_sc['T2_C'])),
        results_sc['T2_C'],
        time_liquid_out if time_liquid_out is not None else range(len(results_liquid['T2_C'])),
        results_liquid['T2_C'],
        time_sc_wht, T_sc_wht_C, time_liquid_wht, T_liquid_wht_C,
        case_sc_name="SC Case",
        case_liquid_name="Liquid Case",
        filename="JT_outlet_temperature_comparison.png"
    )
    
    # --------------------------------------------------
    # FINAL SUMMARY
    # --------------------------------------------------
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nSummary:")
    print(f"  SC CO2 outlet temp range: {np.nanmin(results_sc['T2_C']):.2f}°C to {np.nanmax(results_sc['T2_C']):.2f}°C")
    print(f"  Liquid CO2 outlet temp range: {np.nanmin(results_liquid['T2_C']):.2f}°C to {np.nanmax(results_liquid['T2_C']):.2f}°C")
    
    # Phase analysis
    sc_phases = set([p for p in results_sc['phases'] if p != "FAILED"])
    liquid_phases = set([p for p in results_liquid['phases'] if p != "FAILED"])
    print(f"\nPhases encountered:")
    print(f"  SC case: {sc_phases}")
    print(f"  Liquid case: {liquid_phases}")
    
    # Two-phase analysis
    sc_two_phase = np.sum(~np.isnan(results_sc['Q']))
    liquid_two_phase = np.sum(~np.isnan(results_liquid['Q']))
    print(f"\nTwo-phase points:")
    print(f"  SC case: {sc_two_phase} points")
    print(f"  Liquid case: {liquid_two_phase} points")


if __name__ == "__main__":
    main()
