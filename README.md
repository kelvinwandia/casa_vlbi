# HSA (VLBA + Y1) ANTAB and Calibration Pipeline

This repository provides a Python-based pipeline for processing **HSA (VLBA + Y1)** data using **AIPS**, **ParselTongue**, and **CASA**.  
It automates the process of loading FITS-IDI files into AIPS, generating antenna tables (ANTAB) files, attaching system temperatures and gain curves using scripts from JIVE and performing standard VLBI data calibration steps

---

## 📘 Overview

This pipeline performs the following operations:

1. **Load FITS-IDI data** into AIPS using ParselTongue.  
2. **Generate ANTAB tables** containing system temperature and gain curve information.  
3. **Download and apply JIVE scripts** to attach system and gain curves to the FITS files  
4. **Download and apply JIVE scripts** 


---

## 🧩 Dependencies

### Core Requirements

| Software | Purpose | Installation |
|-----------|----------|---------------|
| **AIPS** | Astronomical Image Processing System (installed system-wide) | [NRAO AIPS](https://www.aips.nrao.edu/) |
| **ParselTongue** | Python interface for AIPS | Installed inside a conda environment | https://www.jive.eu/jivewiki/doku.php?id=parseltongue:parseltongue
| **CASA** | Common Astronomy Software Applications | Installed inside a conda environment |

### Recommended Environment Setup

AIPS should be installed **system-wide**, while **ParselTongue** and **CASA** should be installed in a **conda environment**.
