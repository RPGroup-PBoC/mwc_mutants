---
status: Accepted 
notes: Lowest repressor copy number is above 1 -- interesting
---
# 2019-03-08 Y20I-Q294K IPTG titration

## Purpose
This experiment measures the fold-change in gene expression for the DNA binding 
mutant Y20I across several repressor copy numbers 

## Strain Information

| Location | Plasmid | Genotype | Host Strain | Shorthand | Class |
| :------- | :------ | :------- | ----------: | --------: | -----:| 
| `Manuel 01 - F2` | `pZS4*5-mCherry` |  `galK<>KD4Kan` | `HG104` | `auto` | `auto`|
| `Manuel 02 - F4` | `pZS4*5-mCherry` | `galK<>25-O3+11-YFP ΔlacIZYA` | `HG105` | `delta`| `delta`|

## Titration Series
| Inducer | Concentration |
| :------ | ------------: |
| Isopropylthiogalactopyranoside (IPTG) | 0, 0.1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, 5000  [mM] |

## Analysis Files

**Induction Profile**
![](output/20190308_r1_Y20i_fold_change_curve.png)

## Experimental Protocol

### Cell Husbandry

1. Cells as described in "Strain Information" were grown to saturation overnight in 1mL of LB Miller + spectinomycin in a 2mL-deep 96 well plate.

2. A 96 well plate with 495µL of M9 + 0.5% glucose and appropriate IPTG concentration was prepared via dilution from aliquoted 100X IPTG stocks. 90µL of each  was transferred to a shallow 96 well plate for flow cytometry and kept at 4°C.

3. Cells were diluted 1:10 into 1mL of fresh LB and thoroughly mixed. This was further diluted 1:100 into the M9 IPTG medium.

4. The 96 well plate was placed in the 37°C incubation room and allowed to grow for 8 hours shaking at 225 RPM. This time is sufficient for cells to reach an  OD<sub>600nm</sub>of approximately 0.3

5.  Once 8 hours had passed, the cells were diluted 1:10 in to the chilled shallow 96 well plate filled with media and the appropriate concentration of IPTG. This plate was then taken for measurement via flow cytometry.


**Flow Cytometry**
1. Flow cytometer was cleaned with 1% bleach followed by a flush with washing buffer. If needed, the instrument was calibrated using the MACSQuant calibration beads.
 
2. The diluted 96 well plate was placed on a 4°C ice block for the duration of the measurement. 

2. Aliquots of 80 µL were withdrawn from each well, gently mixed, and analyzed for 100,000
individual events. The voltage settings of the PMTs were as follows:

| Wavelength | Channel | Sensor Voltage|
|:---|:---|---:|
| 488 nm | Forward Scatter (FSC) | 423 V|
| 488 nm   | Side Scatter (SSC) | 537 V|
| 488 nm | Intensity (B1 Filter, 525/50 nm) | 790 V|
| 488 nm | Trigger (debris threshold) | 24.5 V|

3. The measurement was periodically monitored and the buffer volumes were topped
off as necessary.

4. Once completed, the data was collected and transferred to the backup server.