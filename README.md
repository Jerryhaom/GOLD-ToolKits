# GOLD-ToolKits
To promote the applications of Gompertz law-based (GOLD) algorithm of biological age, we developed the R package 'GOLD-ToolKits' as an integration of these approaches.

# Algorithms
| Method | Data Requirements | Primary Function | Key Innovation |
|--------|------------------|------------------|----------------|
|**GOLD**| Longitudinal survival data | Estimate mortality hazard and biological age | Gompertz law application for mortality estimation |
|**GOLD-R**| Cross-sectional data | Estimate biological age residuals | Residual-based approach for cross-sectional studies |
|**GOLD-RL**| Longitudinal survival data | Estimate mortality hazard and age residuals | Combines hazard estimation with residual analysis |
|**GOLD-B**| Compatible with various data types | Update biological age residuals | Bayesian framework for residual refinement |

## 📚 Documentation

### Methodology Formulas
- [View Formulas Documentation](https://raw.githack.com/Jerryhaom/GOLD-ToolKits/main/Formulas.html)

# INSTALL <br>
install.packages("devtools") <br>
devtools::install_github("Jerryhaom/GOLD-ToolKits") <br>
browseVignettes("GOLDToolkits") <br>
[View Online Vignette](https://raw.githack.com/Jerryhaom/GOLD-ToolKits/main/doc/my-vignette.html)  <br>

 
# More information

Contact Information: haombio@gmail.com <br>
## Citation
If you use GOLD-Toolkits in your research, please cite:
1. Hao, M., Zhang, H., et al. (2025). Gompertz Law-Based Biological Age (GOLD BioAge): A Simple and Practical Measurement of Biological Ageing. *Advanced Science*.
2. Zhang, H., Zhang, S., et al. (2025). A Residual Approach to Estimate Biological Age from Gompertz Modeling. *Under Review*.
