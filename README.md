# SRT-Server: Powering the analysis of spatially resolved transcriptomics studies <br>

## Overview
We develop SRT-Server to perform the analysis of spatial resolved transcriptomics (SRT) data. It includes quality control (QC), spatially variable gene detection (SVG), deconvolution (DECON and DECON_PY), cell typing (CT), spatial domain detection (SDD), differential gene expression analysis (DEG), gene set enrichment analysis (GSEA), cell-cell communication identification (CCC) and pseudo-time trajectory inference (TRAJ). <br>
The user can use the server to visit the website: https://spatialtranscriptomicsanalysis.com

## Update log
+ <strong>Version 3.1 (Nov 2023)</strong><br>
We update the DECON_PY and TRAJ module, add the Tutorial webpage and Example Data webpage, and debug some errors. <br>
+ For Example Data webpage, we add the demo video for quick start and three case studies and add the Case Study 2.<br>
+ For Debug buttom, we add the result downloding.<br>
+ <strong>Version 3.0 (Oct 2023)</strong><br>
We update the DECON_PY and TRAJ module, add the Tutorial webpage and Example Data webpage, and debug some errors. <br>
  + For DECON_PY, we add two methods: cell2location and Tangram.<br>
  + For TRAJ, we set DECON as another upstream module, which estimates the pseudo-time trajectory in one cell type. <br>
  + For Tutorial webpage, we included three parts to describe each module: i) a brief introduction and description of the module; ii) the analysis component; iii) the plot component.<br>
  + For Example Data webpage, we included three parts: i) Quick Start; ii) Case Study 1; iii) Case Study 3. In the two case studies, we added Demo Video to show the analytic procedure and parameter setting. <br>
+ <strong>Version 2.2 (Apr 2023)</strong><br>
We update the reference panels of single cell RNA-seq data and marker data.  Totally, we collect 51 datasets for DECON and CT-Annot. We fix some bugs in DECON and CCC_Plot modules. We deploy the server on Tencent Cloud. <br>
+ <strong>Version 2.1 (Feb 2023)</strong><br>
We update R into version 4.2.2. We added the “Without login” function for SRT-Server. SRT-Server in “Without login” setting only give user the download link of result, and the analysis flow chart and parameter settings are not saved when the page closed. <br>
+ <strong>Version 2.0 (Jan 2023)</strong><br>
We update R into version 4.2.1 and fix some bugs. We combine the Annot module to CT module. <br>
+ <strong>Version 1.0 (Dec 2022)</strong><br>
We first build the SRT-Server to power the spatial resolved transcriptomics (SRT) data analysis. 

## Q&A
If there are any problems for SRT-Server, please feel free to ask here. We will ask you as soon as possible. 

## Citation
If you used our server, please cite the our paper: <br>
<em> Sheng Yang, Xiang Zhou, 2023, <strong> SRT-Server: Powering the analysis of spatially resolved transcriptomics studies. </strong> </em>
