# -*- coding: utf-8 -*-
"""
Created on Wed May 21 16:47:44 2025

@author: jose.pauyac
"""



#Tutorial
# https://en.wikipedia.org/wiki/Cubic_equations_of_state
# https://www.youtube.com/watch?v=vi44dy7OajI

# PR EOS

import numpy as np
#from scipy.misc import derivative

#----------------------------------------------------------------------------

#Wellbore Water



#----------------------------------------------------------------------------

# Reservoir Water

#---------------------------------

# T2WELL/ECO2N - IFP67

'''

    Water properties in TOUGH2/ECO2N are calculated from the steam table 
    equations as given by the International Formulation Commitee (1967) IFC-67

    IFC-67 Formulation

    It uses the Tumlirz-Tammann Equation Correlation for water density
    It uses the Kestin Equation Correlation for water viscosity


'''

class IFC67:
    
    def IFC67_water_density(T):
        """
        Approximate liquid water density [kg/m³] from IFC-67 formulation.
        Valid from ~273 K to ~623 K.
        """
        if T < 273.15 or T > 623.15:
            raise ValueError("Temperature out of range for IFC-67 water density (273–623 K).")
        
        # Empirical constants
        a1 = 999.83952
        a2 = 16.945176
        a3 = -7.9870401e-3
        a4 = -46.170461e-6
        a5 = 105.56302e-9
        a6 = -280.54253e-12
        a7 = 0.01687985
    
        t = T - 273.15  # Convert to °C
    
        # IFC-67 density formula (in kg/m³)
        IFC67_rho_water = a1 + a2 * t + a3 * t**2 + a4 * t**3 + a5 * t**4 + a6 * t**5
        IFC67_rho_water /= (1 + a7 * t)
        
        return IFC67_rho_water  # kg/m³
    
    def viscosity_IFC67(T):
        """
        Approximate dynamic viscosity [Pa·s] for water using IFC-67-compatible fit.
        Valid roughly from 273 K to 623 K.
        """
        if T < 273.15 or T > 623.15:
            raise ValueError("Temperature out of range for IFC-67 water viscosity (273–623 K).")
    
        # Viscosity fit from Kestin et al., for liquid water
        A = 2.414e-5  # Pa·s
        B = 247.8     # K
        C = 140       # K
    
        IFC67_mu_water = A * 10**(B / (T - C))  # Pa·s
        return IFC67_mu_water

#---------------------------------

class GEM_correlations_water:
    
    def rowe_chou_water_density(T, P):
        """
        Compute water density [kg/m³] using Rowe & Chou (1970) correlation.
        
        Parameters:
        - T_C: Temperature in °C
        - P_Pa: Pressure in Pa

        Returns:
        - Density [kg/m³]
        """
        # Saturation density at 1 atm [kg/m³] — empirical fit from IAPWS
        rho0 = (999.83952 + 16.945176 * T - 7.9870401e-3 * T**2 
               - 46.170461e-6 * T**3 + 105.56302e-9 * T**4 
               - 280.54253e-12 * T**5) / (1 + 16.87985e-3 * T)

        # Compressibility [1/Pa] — valid for 0–100°C
        beta = (5.1e-10 + 3e-13 * (T - 25))  # 1/Pa

        # Atmospheric pressure in Pa
        P0 = 101325.0

        # Corrected density at pressure P
        GEM_correlation_rho_water = rho0 / (1 - beta * (P - P0)) # kg/m3
        return GEM_correlation_rho_water
    
    def kestin_water_viscosity(T):
        """
        Compute dynamic viscosity of water [Pa·s] using the Kestin et al. (1978) correlation.
        
        Parameters:
        - T_C: Temperature in Celsius (scalar or array)
    
        Returns:
        - Viscosity in Pa·s
        """
        T_K = T + 273.15  # Convert to Kelvin
    
        A = 2.414e-5  # Pa·s
        B = 247.8     # K
        C = 140       # K
    
        GEM_correlation_mu_water = A * 10**(B / (T_K - C))  # Pa·s
    
        return GEM_correlation_mu_water

#---------------------------------

from CoolProp.CoolProp import PropsSI

class NIST_water_properties:
    
    def NIST_water_density(T, P):
        return PropsSI('D', 'T', T, 'P', P, 'Water')  # kg/m3
    
    def NIST_water_viscosity(T, P):
        return PropsSI('V', 'T', T, 'P', P, 'Water')  # Pa*s

#----------------------------------------------------------------------------

# Plotting

# Setting up the PT grid
P_range = np.linspace(0.1e6, 0.6e8, 100)   #Pa

# Temperatures in Celsius and Kelvin
T_liquid_C = [30, 60, 90]
T_liquid_K = [t + 273.15 for t in T_liquid_C]

T_vapor_C = [110, 150, 200]
T_vapor_K = [t + 273.15 for t in T_vapor_C]

# Helper to collect property curves
def get_water_properties(T_K_list, prop_fn, pressure_dependent=True):
    results = []
    for T in T_K_list:
        prop_vs_P = []
        for P in P_range:
            try:
                if pressure_dependent:
                    val = prop_fn(T, P)
                else:
                    val = prop_fn(T)
            except:
                val = np.nan
            prop_vs_P.append(val)
        results.append(prop_vs_P)
    return results


# Define data for liquid
rho_ifc67_liq = get_water_properties(T_liquid_K, IFC67.IFC67_water_density, pressure_dependent=False)
rho_gem_correlation_liq = get_water_properties(T_liquid_C, GEM_correlations_water.rowe_chou_water_density, pressure_dependent=True)
rho_nist_liq = get_water_properties(T_liquid_K, NIST_water_properties.NIST_water_density)

mu_ifc67_liq = get_water_properties(T_liquid_K, IFC67.viscosity_IFC67, pressure_dependent=False)
mu_gem_correlation_liq = get_water_properties(T_liquid_C, GEM_correlations_water.kestin_water_viscosity, pressure_dependent=False)
mu_nist_liq = get_water_properties(T_liquid_K, NIST_water_properties.NIST_water_viscosity)

#---------------------------------

# # Separate plots

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Reservoir CO2

#---------------------------------

# T2WELL/ECO2N - CO2TAB

'''

CO2 properties in TOUGH2/ECO2N for pure CO2 are 
obtained from correlationsdeveloped by Altunin et al. (1975), 
presented in a tabulated from in the CO2TAB file


'''

import pandas as pd

class CO2TAB_EOS:
    supercritical_table = None
    gaseous_table = None
    liquid_table = None

    @classmethod
    def load_excel(cls, file_path):
        df_full = pd.read_excel(file_path, sheet_name=0, header=None)

        df_supercritical = df_full.iloc[5:27, 1:8]
        df_gaseous = df_full.iloc[5:27, 11:18]
        df_liquid = df_full.iloc[5:28, 21:28]

        cls.supercritical_table = cls._process_table(df_supercritical)
        cls.gaseous_table = cls._process_table(df_gaseous)
        cls.liquid_table = cls._process_table(df_liquid)
        
        return df_supercritical, df_gaseous, df_liquid

    @staticmethod
    def _process_table(df_slice):
        df = df_slice.copy()
        df.columns = df.iloc[0]
        df = df[1:].reset_index(drop=True)
        cols = ['T, degC', 'P, bar', 'Density, kg/m3', 'Viscosity, Pa*s', 'Enthalpy, J/kg', 'Viscosity, cP', 'Enthalpy, kJ/kg']
        df[cols] = df[cols].apply(pd.to_numeric, errors='coerce')
        return df
    
    @classmethod
    def get_density(cls, T, P):
        """Get density (kg/m³) from CO2TAB tables at given T (K) and P (Pa)."""
        T_degC = T - 273.15
        P_bar = P / 1e5  # Convert Pa to bar
    
        df = cls._get_phase_table(T_degC, P_bar)
        if df is None:
            return np.nan
        idx = (df['P, bar'] - P_bar).abs().idxmin()
        return df.iloc[idx]['Density, kg/m3']
    
    @classmethod
    def get_viscosity(cls, T, P, unit='Pa*s'):
        T_degC = T - 273.15
        P_bar = P / 1e5
        df = cls._get_phase_table(T_degC, P_bar)
        if df is None:
            return np.nan
        idx = (df['P, bar'] - P_bar).abs().idxmin()
        return df.iloc[idx]['Viscosity, cP'] if unit == 'cP' else df.iloc[idx]['Viscosity, Pa*s']
    
    # @classmethod
    # def get_enthalpy(cls, T, P, unit='J/kg'):
    #     """Get enthalpy (kJ/kg) from CO2TAB tables at given T (K) and P (Pa)."""
    #     T_degC = T - 273.15
    #     P_bar = P / 1e5
    #     df = cls._get_phase_table(T_degC, P_bar)
    #     if df is None:
    #         return np.nan
    #     idx = (df['P, bar'] - P_bar).abs().idxmin()
    #     # Convert J/kg to kJ/kg
    #     return df.iloc[idx]['Enthalpy, kJ/kg'] if unit == 'kJ/kg' else df.iloc[idx]['Enthalpy, J/kg']

    @classmethod
    def _get_phase_table(cls, T_degC, P_bar):
        """Determine which phase table to use based on T and P."""
        if cls.supercritical_table is not None:
            df_sc = cls.supercritical_table
            mask_sc = (df_sc['T, degC'] == T_degC) & (df_sc['P, bar'] >= P_bar)
            if not df_sc[mask_sc].empty:
                return df_sc
        if cls.gaseous_table is not None:
            df_gas = cls.gaseous_table
            mask_gas = (df_gas['T, degC'] == T_degC) & (df_gas['P, bar'] >= P_bar)
            if not df_gas[mask_gas].empty:
                return df_gas
        if cls.liquid_table is not None:
            df_liq = cls.liquid_table
            mask_liq = (df_liq['T, degC'] == T_degC) & (df_liq['P, bar'] >= P_bar)
            if not df_liq[mask_liq].empty:
                return df_liq
        return None

    @classmethod
    def get_property_data(cls, property_name, T_degC, phase=None):
        phases = {
            'Supercritical': cls.supercritical_table,
            'Gaseous': cls.gaseous_table,
            'Liquid': cls.liquid_table
        }
    
        if phase:
            df = phases.get(phase)
            if df is not None:
                filtered = df[df['T, degC'] == T_degC]
                if not filtered.empty:
                    return filtered['P, bar'].values, filtered[property_name].values
        else:
            for df in phases.values():
                filtered = df[df['T, degC'] == T_degC]
                if not filtered.empty:
                    return filtered['P, bar'].values, filtered[property_name].values
    
        return None, None
        
# CO2TAB_EOS.load_excel('T2WELL_EOS_CO2TAB.xlsx')

#---------------------------------

# GEM: PR

class PR_EOS:
    
    @staticmethod
    def PR_CO2(T, P, phase='vapor'):
        """
        Peng-Robinson EOS to calculate compressibility factor Z.

        Args:
            T (float): Temperature in Kelvin.
            P (float): Pressure in Pascals.

        Returns:
            float: Compressibility factor Z.
        """
        
        # Critical properties of CO2
        T_c = 304.12  # K
        P_c = 73.77e5  # Pa
        R = 8.314  # J/(mol·K)
        omega = 0.228 # CO2 acentric factor
        
        # EOS parameters
        a = 0.45724 * R**2 * T_c**2 / P_c
        b = 0.07780 * R * T_c / P_c
        kappa = 0.37464 + 1.54226 * omega - 0.26992 * (omega)**2

        def alpha_CO2(T):
            T_r = T / T_c
            return (1 + kappa * (1 - np.sqrt(T_r))) ** 2
        
        A = a * alpha_CO2(T) * P / (R**2 * T**2)
        B = b * P / (R * T)
        
        # Cubic equation coefficients
        coeffs = [1, -(1 - B), (A - 3 * (B ** 2) - 2 * B), -(A * B - B ** 2 - B ** 3)]
        
        # Solve cubic equation for Z
        roots = np.roots(coeffs)
    
        # Filter real roots
        real_roots = [r.real for r in roots if np.isreal(r)]
        real_roots = sorted([r.real for r in roots if np.isreal(r)])
        if len(real_roots) == 0:
            raise ValueError("No real roots found for the cubic equation.")
        
        if phase == 'vapor':
            PR_Z_CO2 = max(real_roots)  # largest root = vapor
        elif phase == 'liquid':
            PR_Z_CO2 = min(real_roots)  # smallest root = liquid
        else:
            raise ValueError("phase must be 'vapor' or 'liquid'")
        
        return PR_Z_CO2
    
class PR_CO2_properties:
    
    @staticmethod
    def PR_CO2_compressibility_coefficient(T, P):
        
        """
        Calculate the compressibility coefficient of CO2.

        Args:
            T (float): Temperature in Kelvin.
            P (float): Pressure in Pascals.

        Returns:
            float: Compressibility coefficient (1/Pa).
        """
        Z = PR_EOS.PR_CO2(T, P)
        
        R = 8.314  # J/(mol·K)
        # M_CO2 = 44.01e-3  # kg/mol (molar mass of CO2)
        
        # Molar volume (m^3/mol)
        Vm = Z * R * T / P
        
        # Numerical differentiation of Vm with respect to P
        delta_P = 1e-5 * P  # Small perturbation
        Z_perturbed = PR_EOS.PR_CO2(T, P + delta_P)
        Vm_perturbed = Z_perturbed * R * T / (P + delta_P)
        dVdP = (Vm_perturbed - Vm) / delta_P
        
        # Compressibility coefficient (1/Pa)
        c = -dVdP / (Vm+1e-5)
        
        return c
    
    def PR_CO2_density(T, P, phase='vapor'):
        
        Z = PR_EOS.PR_CO2(T, P)
        R = 8.314  # J/(mol·K)
        # M_CO2 = 44.01e-3  # kg/mol (molar mass of CO2)
        
        # Molar volume (m^3/mol)
        Vm = Z * R * T / P
        
        # Molar density
        rho_molar = 1 / (Vm+1e-5) # mol/m3
        
        # Molar mass
        M = 44.009 # g/mol
        
        # Mass density
        PR_rho_CO2 =  rho_molar * (M / 1000) #kg/m^3
        
        return PR_rho_CO2
    
    def herning_zipperer_viscosity_pure_CO2(T):
        
        """
        Calculate viscosity of pure CO2 using Herning and Zipperer structure.
        
        Parameters:
        - T (float or np.ndarray): Temperature in Kelvin
    
        Returns:
        - mu_CO2 (float or np.ndarray): Viscosity of CO2 in Pa·s
        """
        
        # Pure CO2 viscosity model (Sutherland-like or simplified empirical fit)
        # This is an empirical approximation similar to NIST's formulation
        # You can replace this with Vesovic et al. or other detailed models.
        
        # Constants (fitted to CO2 data in gas phase)
        A = 1.37e-6   # [Pa·s]
        B = 240       # [K]
        HZYT_mu_CO2 = A * (T ** 1.5) / (T + B)
    
        return HZYT_mu_CO2 * 1000 # cP
    
    def PR_CO2_enthalpy(T, P, phase='vapor'):
        """
        Calculate the enthalpy of CO2 using the Peng-Robinson EOS.
    
        Args:
            T (float): Temperature in Kelvin.
            P (float): Pressure in Pascals.
    
        Returns:
            float: Enthalpy in kJ/kg
        """
        # Constants
        T_c = 304.12  # K
        P_c = 73.77e5  # Pa
        R = 8.314  # J/(mol·K)
        omega = 0.228  # CO2 acentric factor
        M_CO2 = 44.009e-3  # kg/mol
    
        # EOS parameters
        a = 0.45724 * R**2 * T_c**2 / P_c
        b = 0.07780 * R * T_c / P_c
        kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega ** 2
    
        T_r = T / T_c
        alpha = (1 + kappa * (1 - np.sqrt(T_r))) ** 2
        da_dT = -0.45724 * R**2 * T_c**2 / P_c * (
            kappa * (1 / np.sqrt(T_r)) * (1 / T_c) * (1 + kappa * (1 - np.sqrt(T_r)))
        )
    
        aT = a * alpha
        Z = PR_EOS.PR_CO2(T, P, phase)
        B = b * P / (R * T)
    
        # Departure enthalpy in J/mol
        ln_term = np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B))
        h_dep = R * T * (Z - 1) + ((T * da_dT - aT) / (2 * np.sqrt(2) * b)) * ln_term
    
        # Ideal gas enthalpy approximation (can be replaced by cp integral)
        cp_ideal = 37.135  # J/mol·K, average for CO2 gas in 250-600K
        h_ideal = cp_ideal * (T - 298.15)  # reference at 298.15 K (25C)
    
        # Total molar enthalpy (J/mol)
        h_total_mol = h_ideal + h_dep
    
        # Convert to kJ/kg
        PR_h_CO2 = h_total_mol / M_CO2  # J/kg
    
        return PR_h_CO2 / 1000 # kJ/kg
    
# --- CoolProp-backed PR properties for CO2 ---
class PR_CO2_properties_CoolProp:
    """CO2 properties using CoolProp's Peng–Robinson backend.
       Density/Enthalpy via PR, viscosity via HEOS (PR has no transport)."""

    R_univ = 8.31446261815324       # J/mol/K
    M_CO2  = 44.0095e-3             # kg/mol
    R_spec = R_univ / M_CO2         # J/kg/K

    @staticmethod
    def PR_CO2_density(T, P, phase=None):
        """Density [kg/m3] with PR backend."""
        return PropsSI('D', 'T', T, 'P', P, 'PR::CO2')

    @staticmethod
    def PR_CO2_enthalpy(T, P, phase=None):
        """Enthalpy [kJ/kg] with PR backend."""
        return PropsSI('H', 'T', T, 'P', P, 'PR::CO2') / 1000.0

    @staticmethod
    def herning_zipperer_viscosity_pure_CO2(T):
        
        """
        Calculate viscosity of pure CO2 using Herning and Zipperer structure.
        
        Parameters:
        - T (float or np.ndarray): Temperature in Kelvin
    
        Returns:
        - mu_CO2 (float or np.ndarray): Viscosity of CO2 in Pa·s
        """
        
        # Pure CO2 viscosity model (Sutherland-like or simplified empirical fit)
        # This is an empirical approximation similar to NIST's formulation
        # You can replace this with Vesovic et al. or other detailed models.
        
        # Constants (fitted to CO2 data in gas phase)
        A = 1.37e-6   # [Pa·s]
        B = 240       # [K]
        HZYT_mu_CO2 = A * (T ** 1.5) / (T + B)
    
        return HZYT_mu_CO2 * 1000 # cP

    @staticmethod
    def PR_CO2_Z(T, P, phase=None):
        """Compressibility factor Z from PR density."""
        rho = PR_CO2_properties_CoolProp.PR_CO2_density(T, P)
        return P / (rho * PR_CO2_properties_CoolProp.R_spec * T)

#---------------------------------

# NIST

class NIST_CO2_properties:
    
    def NIST_CO2_density(T, P):
        return PropsSI('D', 'T', T, 'P', P, 'CO2')  # kg/m3
    
    def NIST_CO2_viscosity(T, P):
        return PropsSI('V', 'T', T, 'P', P, 'CO2') * 1000  # cP
    
    def NIST_CO2_enthalpy(T, P):
        return PropsSI('H', 'T', T, 'P', P, 'CO2')  / 1000 # kJ/kg

    def NIST_CO2_saturation_pressure(T):
        T_crit = 304.12
        if T >= T_crit:
            return np.inf  # Supercritical: no saturation pressure
        return PropsSI('P', 'T', T, 'Q', 0, 'CO2')  # Pa

#----------------------------------------------------------------------------

# def get_saturation_pressure_CO2(T_K):
#     """
#     Returns the saturation pressure (Pa) of CO2 at temperature T_K.
#     Returns infinity if T_K is above the critical temperature.
#     """
#     T_critical = 304.12  # Critical temperature of CO2 in K
#     if T_K >= T_critical:
#         return np.inf  # No liquid phase exists above critical temperature
#     try:
#         P_sat = PropsSI('P', 'T', T_K, 'Q', 0, 'CO2')  # Saturation pressure in Pa
#         return P_sat
#     except:
#         return np.nan

#----------------------------------------------------------------------------

def get_saturation_pressure_CO2(T_K):
    try:
        return NIST_CO2_properties.NIST_CO2_saturation_pressure(T_K)  # Returns Pa
    except:
        return np.nan

def get_gaseous_CO2_properties(T_K_list, prop_fn, pressure_dependent=True):
    values_per_T = []
    pressures_per_T = []

    for T in T_K_list:
        P_sat = get_saturation_pressure_CO2(T)
        prop_vs_P = []
        valid_P = []

        for P in P_gaseous_CO2:
            if np.isnan(P_sat) or P < P_sat:  # Filter out supercritical/liquid
                try:
                    val = prop_fn(T, P) if pressure_dependent else prop_fn(T)
                except:
                    val = np.nan
                prop_vs_P.append(val)
                valid_P.append(P)

        values_per_T.append(prop_vs_P)
        pressures_per_T.append(valid_P)

    return values_per_T, pressures_per_T

#----------------------------------------------------------------------------

# Plotting

# Setting up the PT grid
P_SC_CO2 = np.linspace(0.76e7, 0.200e8, 100)   #Pa
T_C_SC_CO2 = [41, 81, 101]
T_K_SC_CO2 = [t + 273.15 for t in T_C_SC_CO2]

P_gaseous_CO2 = np.linspace(0.1e6, 0.706e7, 100)   #Pa
T_C_gaseous_CO2 = [5.04, 21.04, 41]
T_K_gaseous_CO2 = [t + 273.15 for t in T_C_gaseous_CO2]

# Setting up the PT grid
#P_liquid_CO2 = np.linspace(0.48536e7, 0.476e8, 100)
P_liquid_CO2 = np.linspace(0.585e7, 0.200e8, 100)   #Pa
T_C_liquid_CO2 = [11.04, 21.04, 29.04]
T_K_liquid_CO2 = [t + 273.15 for t in T_C_liquid_CO2]

CO2TAB_EOS.load_excel("T2WELL_EOS_CO2TAB_short.xlsx")

# Define linestyles per EOS
linestyles = {
    'PR': '-',
    'NIST': '--',
    'T2WELL': ':'
}

# Define a color map for unique temperatures
unique_T_C_SC_CO2 = T_C_SC_CO2
#colors = cm.viridis(np.linspace(0, 1, len(unique_T_C)))  # Or use 'tab10', 'plasma', etc.
colors = ['blue', 'red', 'green']

# Create a mapping from temperature to color
temp_color_map_SC_CO2 = dict(zip(unique_T_C_SC_CO2, colors))

# Define a color map for unique temperatures
unique_T_C_gaseous_CO2 = T_C_gaseous_CO2

# Create a mapping from temperature to color
temp_color_map_gaseous_CO2 = dict(zip(unique_T_C_gaseous_CO2, colors))

# Define a color map for unique temperatures
unique_T_C_liquid_CO2 = T_C_liquid_CO2

# Create a mapping from temperature to color
temp_color_map_liquid_CO2 = dict(zip(unique_T_C_liquid_CO2, colors))

# Define linestyles per EOS
linestyles_T2WELL_SW = {
    'NIST': '--',
    'T2WELL': ':'
}

# Define a color map for unique temperatures
unique_T_C_SC_CO2_T2WELL_SW = T_C_SC_CO2
#colors = cm.viridis(np.linspace(0, 1, len(unique_T_C)))  # Or use 'tab10', 'plasma', etc.
colors_T2WELL_SW = ['blue', 'red', 'green']

# Create a mapping from temperature to color
temp_color_map_SC_CO2_T2WELL_SW = dict(zip(unique_T_C_SC_CO2_T2WELL_SW, colors_T2WELL_SW))

# Define a color map for unique temperatures
unique_T_C_gaseous_CO2_T2WELL_SW = T_C_gaseous_CO2

# Create a mapping from temperature to color
temp_color_map_gaseous_CO2_T2WELL_SW = dict(zip(unique_T_C_gaseous_CO2_T2WELL_SW, colors_T2WELL_SW))

# Define a color map for unique temperatures
unique_T_C_liquid_CO2_T2WELL_SW = T_C_liquid_CO2

# Create a mapping from temperature to color
temp_color_map_liquid_CO2_T2WELL_SW = dict(zip(unique_T_C_liquid_CO2_T2WELL_SW, colors_T2WELL_SW))

import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(7, 4.5), sharex='col')
plt.subplots_adjust(hspace=0.15, wspace=0.35, bottom=0.21)  # leave room for legends at the top

def plot_property_NIST_T2WELL(ax, title, x_label, y_label,
                              T_K_list, T_C_list, P_list,
                              prop_func_NIST, t2well_property, phase, temp_color_map):
    phase_map = {'gaseous': 'vapor', 'liquid': 'liquid', 'supercritical': 'vapor'}
    phase_pr = phase_map.get(phase.lower())

    for T_K, T_C in zip(T_K_list, T_C_list):
        color = temp_color_map[T_C]
        P_sat = get_saturation_pressure_CO2(T_K)

        if phase.lower() == 'gaseous':
            pressures_filtered = [P for P in P_list if np.isnan(P_sat) or P < P_sat]
        elif phase.lower() == 'liquid':
            pressures_filtered = [P for P in P_list if np.isnan(P_sat) or P >= P_sat]
        else:
            pressures_filtered = P_list  # supercritical

        pressures_bar = np.array(pressures_filtered) / 1e5

        # NIST
        nist_curve = [prop_func_NIST(T_K, P) for P in pressures_filtered]
        ax.plot(pressures_bar, nist_curve,
                linestyle=linestyles_T2WELL_SW['NIST'],
                color=color, label=f'CoolProp {T_C}°C')

        # T2WELL
        P_t2well_bar, prop_t2well = CO2TAB_EOS.get_property_data(t2well_property, T_C, phase=phase)
        if P_t2well_bar is not None:
            ax.plot(P_t2well_bar, prop_t2well,
                    linestyle=linestyles_T2WELL_SW['T2WELL'],
                    color=color, label=f'T2WELL {T_C}°C')

    ax.set_title(title, fontsize=11)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(True)
    # NOTE: no ax.legend() here — we'll add one legend per column below


# Left column (gaseous)
plot_property_NIST_T2WELL(
    axes[0, 0], 'Gaseous', '', 'Density (kg/m³)',
    T_K_gaseous_CO2, T_C_gaseous_CO2, P_gaseous_CO2,
    NIST_CO2_properties.NIST_CO2_density,
    'Density, kg/m3', 'Gaseous', temp_color_map_gaseous_CO2
)
plot_property_NIST_T2WELL(
    axes[1, 0], '', 'Pressure (bar)', 'Enthalpy (kJ/kg)',
    T_K_gaseous_CO2, T_C_gaseous_CO2, P_gaseous_CO2,
    NIST_CO2_properties.NIST_CO2_enthalpy,
    'Enthalpy, kJ/kg', 'Gaseous', temp_color_map_gaseous_CO2
)

# Right column (supercritical)
plot_property_NIST_T2WELL(
    axes[0, 1], 'Supercritical', '', 'Density (kg/m³)',
    T_K_SC_CO2, T_C_SC_CO2, P_SC_CO2,
    NIST_CO2_properties.NIST_CO2_density,
    'Density, kg/m3', 'Supercritical', temp_color_map_SC_CO2
)
plot_property_NIST_T2WELL(
    axes[1, 1], '', 'Pressure (bar)', 'Enthalpy (kJ/kg)',
    T_K_SC_CO2, T_C_SC_CO2, P_SC_CO2,
    NIST_CO2_properties.NIST_CO2_enthalpy,
    'Enthalpy, kJ/kg', 'Supercritical', temp_color_map_SC_CO2
)

# ---------- One legend per column ----------
def col_legend(col_idx, x_anchor):
    handles, labels = [], []
    for ax in axes[:, col_idx].ravel():
        h, l = ax.get_legend_handles_labels()
        for hi, li in zip(h, l):
            if li not in labels:  # dedupe by label text
                handles.append(hi); labels.append(li)
    #ncols = min(len(labels), ncol_max)
    leg = fig.legend(handles, labels,
                     loc='lower center',
                     bbox_to_anchor=(x_anchor, 0.015),  # x across the top (0..1 in figure coords)
                     ncol=3, fontsize=5.4, frameon=True)
    fig.add_artist(leg)

# Place legends centered above each column (x≈0.25 for left, x≈0.75 for right)
col_legend(0, x_anchor=0.25)
col_legend(1, x_anchor=0.7)

# Save and show
plt.savefig('CO2_properties_comparison_CoolProp_NIST_T2WELL.png', dpi=300, bbox_inches='tight')
plt.savefig('CO2_properties_comparison_CoolProp_NIST_T2WELL.eps', format='eps', bbox_inches='tight')
plt.show()

#----------------------------------

fig, axes = plt.subplots(2, 1, figsize=(3, 4))
#plt.subplots_adjust(hspace=0.15, wspace=0.35, bottom=0.21)  # leave room for legends at the top

def plot_property_NIST_T2WELL(ax, title, x_label, y_label,
                              T_K_list, T_C_list, P_list,
                              prop_func_NIST, t2well_property, phase, temp_color_map):
    phase_map = {'gaseous': 'vapor', 'liquid': 'liquid', 'supercritical': 'vapor'}
    phase_pr = phase_map.get(phase.lower())

    for T_K, T_C in zip(T_K_list, T_C_list):
        color = temp_color_map[T_C]
        P_sat = get_saturation_pressure_CO2(T_K)

        if phase.lower() == 'gaseous':
            pressures_filtered = [P for P in P_list if np.isnan(P_sat) or P < P_sat]
        elif phase.lower() == 'liquid':
            pressures_filtered = [P for P in P_list if np.isnan(P_sat) or P >= P_sat]
        else:
            pressures_filtered = P_list  # supercritical

        pressures_bar = np.array(pressures_filtered) / 1e5

        # NIST
        nist_curve = [prop_func_NIST(T_K, P) for P in pressures_filtered]
        ax.plot(pressures_bar, nist_curve,
                linestyle=linestyles_T2WELL_SW['NIST'],
                color=color, label=f'CoolProp {T_C}°C')

        # T2WELL
        P_t2well_bar, prop_t2well = CO2TAB_EOS.get_property_data(t2well_property, T_C, phase=phase)
        if P_t2well_bar is not None:
            ax.plot(P_t2well_bar, prop_t2well,
                    linestyle=linestyles_T2WELL_SW['T2WELL'],
                    color=color, label=f'T2WELL {T_C}°C')

    ax.set_title(title, fontsize=11)
    ax.set_xlabel(x_label, fontsize=10)
    ax.set_ylabel(y_label, fontsize=10)
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(True)
    ax.legend(fontsize=4.5, ncol=3, loc='best')

plot_property_NIST_T2WELL(
   axes[0], '', 'Pressure [bar]', 'Viscosity (cP)',
   T_K_gaseous_CO2, T_C_gaseous_CO2, P_gaseous_CO2,
   NIST_CO2_properties.NIST_CO2_viscosity,
   'Viscosity, cP', 'Gaseous', temp_color_map_gaseous_CO2
)

plot_property_NIST_T2WELL(
   axes[1], '', 'Pressure [bar]', 'Viscosity (cP)',
   T_K_SC_CO2, T_C_SC_CO2, P_SC_CO2,
   NIST_CO2_properties.NIST_CO2_viscosity,
   'Viscosity, cP', 'Supercritical', temp_color_map_SC_CO2
)

# Save and show
plt.savefig('CO2_properties_comparison_CoolProp_NIST_T2WELL_viscosity.png', dpi=300, bbox_inches='tight')
plt.savefig('CO2_properties_comparison_CoolProp_NIST_T2WELL_viscosity.eps', format='eps', bbox_inches='tight')
plt.show()

#----------------------------------------------------------------------------

# For injection temperature only

data_Tinj_T2WELL_density = {
    'Pressure (bar)': [1, 21, 50.9, 67.5, 70.6, 73.9, 76, 80, 84, 88, 92, 96, 100],
    'Density (kg/m3)': [1.73E+00, 4.01E+01, 1.21E+02, 2.01E+02, 2.26E+02, 2.60E+02, 2.90E+02, 4.23E+02, 5.95E+02, 646, 675, 696, 7.13E+02]
}


# Calculate the density for CO2 using Span and Wagner EOS at 35 C and the given pressures
pressures_CO2LINK = np.linspace(1, 100, 100)
densities_CO2LINK = [PropsSI('D', 'P', p*1e5, 'T', 35 + 273.15, 'CO2') for p in pressures_CO2LINK]

plt.figure(figsize=(3, 2))
plt.plot(pressures_CO2LINK, densities_CO2LINK, color='blue', label="CO2LINK @ Tinj = 35 °C")
plt.plot(data_Tinj_T2WELL_density['Pressure (bar)'], data_Tinj_T2WELL_density['Density (kg/m3)'], color='red', label="T2WELL @ Tinj = 35 °C")
plt.ylabel(r'Density (kg/m$3$)', fontsize=7)
plt.xlabel('Pressure (bar)', fontsize=7)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.grid(True)
plt.legend(fontsize=5)
plt.savefig('CO2_properties_comparison_CoolProp_NIST_T2WELL_Wellhead_Conditions.png', dpi=300, bbox_inches='tight')
plt.savefig('CO2_properties_comparison_CoolProp_NIST_T2WELL_Wellhead_Conditions.eps', format='eps', bbox_inches='tight')

plt.show()

#----------------------------------------------------------------------------
