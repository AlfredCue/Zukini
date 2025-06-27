from pathlib import Path
from src.plot_densities_robustness import ZDPlasKinAnalyzer

#This basic script is just to use for verifications

an = ZDPlasKinAnalyzer()        
an.plot_single_species(
    "N^+",
    Path("Data/QtFiles/Qt_2"),
    Path("Data/QtFiles/Qt_2")
)
