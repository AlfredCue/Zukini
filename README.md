# Zukini
Zukini is a Python framework that integrates ZDPlasKin and PumpKin to reduce and validate NOx plasma kinetics. It supports duty-cycle sweeps, robustness checks over gas, temperature, and power, and turns ZDPlasKin into a virtual bench to compare simulations with real experimental data.


# Zukini-1.0: NOx Chemistry Reduction Framework with ZDPlasKin

**Zukini** is a Python-based orchestration framework for analyzing and reducing NOx chemistry using [ZDPlasKin](http://www.zdplaskin.laplace.univ-tlse.fr/index.html). It supports:
- Reaction set optimization
- Parametric simulation runs            (Needs post processing for conclusions)
- Robustness checking across conditions (On revision)

# NOTE
For the moment the optimization algorithm is configures to work with NOx kinemetics file. Still, the series of functions 
serve as a great base to improve the robustness of the tool allowing it to optimize any type of kinetic mechanism.

## 📁 Project Structure

```
Zukini/
├── bin/
│    └──run_optizd.py                           #Script for running the optimization 
│    └──run_parametric.py                       #Script for doing the parametric run(For now it just varies the power input)
│    └──run_robustness_checker.py               #Script to do the robustness check of any new kinet file produced
├──data/
│    └──QtFiles/
│            └──Qt_1/                           #Folder location for reference QT files
│            └──Qt_2/                           #Folder location for reduced simulation QT files
│    └──Reduced_Sets/                           #Folder location for reduced sets and original kinet file
│            └──kinet_original.inp              #Original kinet file
│            └──kinet_filtered_1.inp            #For the robustness check you have to keep these names
│            └──kinet_filtered_2.inp            
│    └──temp/
│            └──power_info/                     #Helper folder used for parametric run
│            └──qt_backup/                      #Helper folder used for Optimization run
├──images/                                      #Any plot created in the processes will be saved here
│    └──comparisons/                            #Plot for the robustness check
├──input/
│    └──global_input.txt                        #Input file (user configurable) for the optimization
├──logs/                                        #Location for all the log outputs of the scripts
├──output/                                      #Location for all outputs (user interest) of the scripts
│    └──experiments/                            #Folder location for the parametric runs results
│    └──robustness_analysis_results/
│    └──kinet_final.inp                         #When optimization runned, the final kinet will be here
├──resources/
│    └──pumpkin-1.4/                            #Location for all Pumpkin files to use it
│    └──ZDPlaskin/                              #Location for all ZDPlaskin files to use it
├──src/                                         #Here is where all the scripts of the tool are and some helper files
│    └──helper_bat_scripts/
│    └──helper_files/
│    └──matlab_power_generator/                 #Matlab .m and simulink models to run the power generator                    
├──README.md
├──requirements.txt

```

---

## 🚀 Setup Instructions

### 1. Install Dependencies

Install required Python libraries using:

```bash
pip install -r requirements.txt
```

Ensure [MATLAB](https://www.mathworks.com/products/matlab.html) is installed and available via command line (batch mode is used).

If using the [MATLAB] option to generate the power input, ensure you have PLECS installed and the conection modules with matlab
Alternatively you could use pre-made data_in.dat files, but the script for the usage of these is yet to be created

Ensure you have one of these compiler options for compiling fortran code: [Gfortran](https://fortran-lang.org/learn/os_setup/install_gfortran/) or [IFort](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html) Second one requires a licence. Make sure that when installing the compiler, it has an entry in the PATHS of the system.

Its highly recommended to use a virtual enviroment to run the project. You can use conda for this(https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) Follow the page for instructions in creating a Conda env

---

## ZDPLASKIN (Input / output explanation)

ZDPlasKin is configured in the master fortran code (dbd0D.F90) here you can modify
the initial concentration of the gases and the values for preassure and temperature. One can modify the constant temperature
to be a value recieved by an external .dat file, like the case of the power input, basically giving a time variant
temperature to further improve the simulation.

As it is now (Version 1.0) the core files for a ZDPlasKin run are:
```
├──> data_in.dat           #Has the values for time, voltage and current
├──> dbd0D.F90             #Is the master code for running ZDPlasKin
├──> kinet.inp             #The file that contains the kinetic mechanism of the simulation
```

For the outputs, in the master code we are calling a subroutine of ZDPlasKin that generates a series of files
with the suffix "qt" that encapsulate all the information of the simulation.

```
├──> qt_conditions_list.txt           #Names of the conditions considered in the simulation
├──> qt_conditions.txt                #Values per time-step of the conditions of the simulation
├──> qt_densities.txt                 #Values per time-step of the density of each specie given in cm^-3
├──> qt_matrix.txt                    #stecheometric matrix that relates the species with the reactions (-1 reactant, 1 product)
├──> qt_rates.txt                     #Values for the reaction rate of each chemical reaction each time-step
├──> qt_reactions_list.txt            #List of all reactions considered/simulated
├──> qt_species_list.txt              #List of all species considered/simulated
```

---

## Single runs (Experiment)

If you want to do a single run experiment, you can run the script (single_run_module.py) by default it run the 
pre-processor, compiles the code, and runs the simulation with the files located in Zukini/resources/ZDPlaskin. The
output is configured to go to Zukini/data/Qt_files/Qt_2.

### Main Entry:
```bash
python single_run_module.py
```



---


## 🔧 Block 1: Kinetic File Optimization

### Purpose:
Optimize a chemical mechanism by filtering reactions based on:
- Pathway relevance (via Pumpkin and Sensitibity)
Sensitibity could be evaded and still get good results but this measurment could be used for further scientific anlysis
- Reaction rates (group thresholds) This is more for numerical stability

### Main Entry:
```bash
python -m bin.run_optizd
```

### Workflow:
1. Clean and prepare kinetic input. (If you have important qtfiles you want to keep from a single run be sure to save before running this script)
2. Run reference simulation with full model. 
3. Build reactions tree from species of interest with:
   - `sensitivity_filter.py`
   - `pumpkin`
4. Iteratively reduce, test, and recover reactions using:
   - `Filter_1_threshold.py`
   - `brute_force_recovery.py`
   - `Error_calc.py`
5. Save `kinet-final.inp` with a minimal trusted set.

---

## ⚙️ Block 2: Parametric Runs

### Purpose:
Automate simulations across different **duty cycles** to explore system behavior.

### Main Entry:
```bash
python -m bin.run_parametric
```

### Workflow:
- For each duty cycle:
  1. Generate power profile using MATLAB.
  2. Run ZDPlasKin simulation.
  3. Save Qt output and analyze de-NOx performance with `DensityAnalyzer`.

### Output:
Structured experiment folders under `/Experiments/`, including:
- `run_info.txt`
- `qt_*` files
- `period_deltas.csv`, `summary.txt`, and species plots.

---

## 🧪 Block 3: Robustness Checker

### Purpose:
Evaluate how the reduced mechanism performs across varying:
- Temperatures
- Gas mixtures (NO, O₂)
- Duty cycles

### Main Entry:
```bash
python -m bin.run_robustness_checker
```

### Workflow:
1. Loop over all T × Gas × Duty configs.
2. Update gas + temp in Fortran source.
3. Generate power input.
4. Run full and reduced models (original + filtered1 + filtered2).
5. Compare outputs via:
   - `plot_densities_robustness.py`
   - Trust score (R² + NRMSE)
   - Δ Efficiency
6. Save visual + text summaries.

### Output:
Structured under `/robustness_analysis_results/`, includes:
- PNG comparisons
- `metrics_summary.txt`
- Global trust CSV: `compiled_results.csv`

---

## 🧩 Dependencies

- Python 3.9+
- MATLAB (for power profile generation)
- ZDPlasKin compiled via gfortran or ifort
- Compiler for Fortran code, either Gfortran of Ifort
- PLECS for usage in Matlab script
- See `requirements.txt` for Python dependencies

---

## 🧠 Credits

Developed by Alfredo Cuellar Valencia during internship in LAPLACE Laboratory under the tutoring of Professor Hubert Piquet and Doctor Nicolas Bente

---

## 📬 Questions or Extensions?

Feel free to reach out or extend the pipeline with:
- Additional filtering strategies
- GUI integration
- Statistical uncertainty propagation
- More robust usage and runs for ZDPlaskin

## Future works
- Modify the fortran main file so it can accept temperature as an input file, letting ZDPlaskin run simulation with temperature profiles rather than constants
- Build a GUI for easier usage of the tool
-Incorporate QTPlaskin software as the visualization tool

## Contact information

Alfredo Cuellar Valencia
alcuellar4120@gmail.com
