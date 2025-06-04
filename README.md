# Groundwater Flow Modeling By Excel
Excel sheets to solve steady-state groundwater flow for a variety of cases.
The solution of the steady-state groundwater flow equation has been implemented using Excel sheets for the following conceptual models:
- Flow in the horizontal plane in a confined aquifer of known transmissivities
- Flow in the horizontal plane in an unconfined aquifer of known hydraulic conductivities
- Flow in a vertical cross-section of an unconfined aquifer with calculation of the phreatic surface
- In 2025, a new solution to handle transient flow in an unconfined aquifer has been developed

The first case is described in the paper by Gómez-Hernandez (2022),  the second and third ones in the paper by Gómez-Hernández and Secci (2023) and the transient one in the paper by Gómez-Hernández and Secci (2025)

The subfolders with flopy scripts contain Python scripts and input files to replicate the Excel model using MODFLOW with the flopy library on the Python environment. They can be used to verify the accuracy of the spreadsheet results

# References
- J. Jaime Gómez-Hernández, [Teaching Numerical Groundwater Flow Modeling with Spreadsheets](https://doi.org/10.1007/s11004-022-10002-4), _Mathematical Geosciences, 54_, 1121-1138, [doi:10.1007/s11004-022-10002-4](http://doi.org//10.1007/s11004-022-10002-4), 2022.
- J. Jaime Gómez-Hernández and Daniele Secci, [Teaching Numerical Groundwater Flow Modeling with Spreadsheets: Unconfined Aquifers and Vertical Cross-Sections](https://jgomez.webs.upv.es/wordpress/wp-content/uploads/2023/10/TeachingNumericalGWFlowWithExcel2_rev2.pdf), _Mathematical Geosciences_, [doi:10.1007/s11004-023-10112-7](https://doi.org/10.1007/s11004-023-10112-7), 2023.
- J. Jaime Gómez-Hernández and Daniele Secci, [Teaching Transient Numerical Groundwater Flow Modeling with Spreadsheet and Large Language Models](https://jgomez.webs.upv.es/wordpress/wp-content/uploads/2025/06/TeachingWithExcelTransientFlowAccepted.pdf), _Mathematical Geosciences_, [doi:](https://doi.org/), accepted, 2025.
