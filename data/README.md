# Data Description

## Framinham Heart Study
```
data("Cholesterol", package="qrLMM")

```

### Format

newid: a numeric vector indicating the subject on which the measurement was made. It represents the subject number in the sample.
ID: a numeric vector indicating the subject on which the measurement was made. It represents the subject number in the population.
cholst: cholesterol level for patient newid.
sex: a dichotomous gender (0=female, 1=male).
age: age of the patient in years.
year: years elapsed since the start of the study to the current measurement.

Source
Zhang, D., & Davidian, M. (2001). Linear mixed models with flexible distributions of random effects for longitudinal data. Biometrics, 57(3), 795-802.

References
https://www.framinghamheartstudy.org/about-fhs/background.php



## Leukemia Survival
```
data("LeukSurv", package="spBayesSurv")

```
### Format

time:	survival time in days
cens:	right censoring status 0=censored, 1=dead
xcoord:	coordinates in x-axis of residence
ycoord:	coordinates in y-axis of residence
age:	age in years
sex:	male=1 female=0
wbc:	white blood cell count at diagnosis, truncated at 500
tpi:	the Townsend score for which higher values indicates less affluent areas
district:	administrative district of residence

Source
Henderson, R., Shimakura, S., and Gorst, D. (2002), Modeling spatial variation in leukemia survival data, Journal of the American Statistical Association, 97(460), 965-972.


## Veteran Lung Cancer Survival

```
data(cancer, package="survival")
```
### Format

inst:	Institution code
time:	Survival time in days
status:	censoring status 1=censored, 2=dead
age:	Age in years
sex:	Male=1 Female=2
ph.ecog:	ECOG performance score as rated by the physician. 0=asymptomatic, 1= symptomatic but completely ambulatory, 2= in bed <50% of the day, 3= in bed > 50% of the day but not bedbound, 4 = bedbound
ph.karno:	Karnofsky performance score (bad=0-good=100) rated by physician
pat.karno:	Karnofsky performance score as rated by patient
meal.cal:	Calories consumed at meals
wt.loss:	Weight loss in last six months (pounds)
Note
The use of 1/2 for alive/dead instead of the usual 0/1 is a historical footnote. For data contained on punch cards, IBM 360 Fortran treated blank as a zero, which led to a policy within the section of Biostatistics to never use "0" as a data value since one could not distinguish it from a missing value. The policy became a habit, as is often the case; and the 1/2 coding endured long beyond the demise of punch cards and Fortran.

Source
Terry Therneau

References
Loprinzi CL. Laurie JA. Wieand HS. Krook JE. Novotny PJ. Kugler JW. Bartel J. Law M. Bateman M. Klatt NE. et al. Prospective evaluation of prognostic variables from patient-completed questionnaires. North Central Cancer Treatment Group. Journal of Clinical Oncology. 12(3):601-7, 1994.

