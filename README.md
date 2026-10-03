# CO₂ Injection Startup Transients

This repository contains Python post-processing scripts supporting a study of CO₂ injection startup transients using CO2LINK, a dynamic controller coupling LedaFlow (wellbore) and CMG-GEM (reservoir).

The analysis includes:

- Validation against T2WELL/ECO2N predictions.
- Thermophysical-property comparisons using CoolProp and T2WELL/ECO2N tables.
- Field-representative liquid and supercritical CO₂ injection cases under mass-rate and pressure control.

> This code supports a manuscript currently under review.
> If you use this work in research, please acknowledge the repository and cite the associated paper when available.

---

## Main Features

- **Thermophysical-property comparison**  
  Compare CO₂ density, enthalpy, and viscosity calculated using CoolProp with values from T2WELL/ECO2N property tables.  
  Location: `Validation/EOS_Comparison/`

- **Transient validation**  
  Compare wellhead pressure and density, bottomhole pressure and temperature, and wellbore velocity between CO2LINK and T2WELL/ECO2N.  
  Location: `Validation/Trend/`

- **Final wellbore profiles**  
  Compare final thermodynamic profiles, including temperature, common-reference specific enthalpy change, specific heat capacity, Joule–Thomson coefficient, density, and viscosity.  
  Location: `Validation/Profile/`

- **Mass-rate-controlled injection**  
  Compare the startup response of liquid and supercritical CO₂ injection cases using wellhead and bottomhole pressure, temperature, density, mass flow rate, and viscosity.  
  Scripts: `TX_Real_field/compare_trend_data_SI_mass_rate.py`

- **Pressure-controlled injection**  
  Examine the corresponding transient response when wellhead pressure is the controlled variable.  
  Script: `TX_Real_field/compare_trend_data_SI_pressure.py`

- **Joule–Thomson valve analysis**  
  Calculate valve-outlet temperature under an isenthalpic expansion assumption and compare it with simulated wellhead temperature.  
  Scripts: `TX_Real_field/JT_coefficient_comparison_SI_*.py`

- **CO₂ phase evolution**  
  Visualize the evolution of pure-CO₂ phase states along the wellbore for both injection cases and control modes. These plots describe pure-CO₂ states, not CO₂–brine phase distributions.  
  Scripts: `TX_Real_field/PT_profile_*.py`

---

## Dependencies

The scripts use Python and the following packages:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `CoolProp`
- `openpyxl`
- `PyTOUGH` — provides `t2listing` for reading T2WELL output.

Install:

```bash
pip install numpy pandas matplotlib scipy CoolProp openpyxl PyTOUGH
```

## Running the Scripts

These scripts process existing simulation outputs; they do not execute the coupled LedaFlow–CMG-GEM simulations.

Before running:

1. Supply the required CSV, TSV, Excel, and T2WELL output files.
2. Update input paths to match your local folders. Some scripts contain machine-specific paths.
3. Set the working directory for scripts using relative filenames.
4. Use separate output directories for different cases to avoid overwriting files.

Example:

```bash
python compare_trend_data_SI_mass_rate.py
```

Figures are exported as PNG and EPS files. Selected scripts also export numerical results as CSV files.

Field temperature inputs are expected in Celsius. The Joule–Thomson scripts explicitly convert their specified Fahrenheit inlet temperatures to Celsius, and thermophysical calculations use Kelvin where required.

For consistent sizing in multi-panel manuscript figures, generate all panels with the same Matplotlib figure-size settings.
