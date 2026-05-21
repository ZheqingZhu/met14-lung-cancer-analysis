# Frontline MET-TKIs outperform systemic therapy in *MET*ex14 NSCLC: a propensity-weighted synthetic control study

This repository contains the dataset, statistical source code, and digitized extraction coordinates for the associated manuscript: "Frontline MET-TKIs outperform systemic therapy in *MET*ex14 NSCLC: a propensity-weighted synthetic control study."

The analytical pipeline is open-sourced to facilitate methodological transparency and independent validation.

## Repository Structure

* **Data/**: Contains the reconstructed patient-level data and cohort baseline matrices.
  * `final_combined_analysis_data_fixed.csv`: The pooled pseudo-individual patient data (pseudo-IPD) including survival times, censoring status, and treatment group assignments.
  * `master_baseline_table_fixed.csv`: Aggregate baseline characteristics used for propensity score simulation and covariate balancing.

* **Code/**: 
  * `main_statistical_analysis.R`: The core statistical pipeline. It executes Monte Carlo baseline simulations, propensity score generation, inverse probability of treatment weighting (IPTW), overlap weighting (OW), restricted mean survival time (RMST) analysis, and E-value sensitivity computation.

* **KM-PoPiGo_Coordinates/**: Contains digitized coordinate data stored in standardized `.json` format for independent extraction validation.

## Reproducibility and KM-PoPiGo Integration

To facilitate independent validation, this repository includes the raw `.json` coordinate files used during the data digitization process. These files are compatible with [KM-PoPiGo](https://kmpopigo.github.io/), a survival data reconstruction web application.

Researchers can reproduce the extraction process and validate the reconstructed IPD through the following steps:

1. Access the [KM-PoPiGo Web Interface](https://kmpopigo.github.io/).
2. Select the **"Load Project"** function and upload any `.json` file provided in the `KM-PoPiGo_Coordinates/` directory.
3. The application will restore the original digitization markers, axis constraints, and number-at-risk intervals.
4. Execute the extraction algorithm to reconstruct the pseudo-IPD and assess the data fidelity.

## Statistical Pipeline

The provided R script (`main_statistical_analysis.R`) implements the following analytical steps:

1. **Covariate Balancing:** Application of average treatment effect on the treated (ATT) weighting for standard systemic therapies, and overlap weighting (OW/ATO) for the immunotherapy cohort (KEYNOTE-024) to address positivity violations.
2. **Survival Analysis:** Cox proportional hazards regression models.
3. **Temporal Dynamics:** Restricted mean survival time (RMST) analysis over a 36-month window to evaluate time-varying treatment effects in the context of non-proportional hazards.
4. **Sensitivity Analysis:** E-value computation to estimate the minimum strength of unmeasured confounding required to attenuate the observed survival effects.

### Requirements

Execution of the R script requires the following packages:
```R
install.packages(c("dplyr", "readr", "survival", "cobalt", "EValue"))
```
## Citation

If you use the data, code, or the JSON constraints from this repository, please cite the original manuscript:
> *(Citation details will be updated upon publication)*

## License

This project is licensed under the [MIT License](LICENSE).
