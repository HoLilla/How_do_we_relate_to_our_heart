# How do we relate to our heart? Neurobehavioral differences across three types of engagement with cardiac interoception

Biological Psychology,
Volume 165,
2021,
108198,
ISSN 0301-0511,
[https://doi.org/10.1016/j.biopsycho.2021.108198.](https://www.sciencedirect.com/science/article/pii/S0301051121001915)

## Abstract: 
Standard measures of interoception are typically limited to the conscious perception of heartbeats. 
However, the fundamental purpose of interoceptive signaling, is to regulate the body. 
We present a novel biofeedback paradigm to explore the neurobehavioral consequences of three different types of 
engagement with cardiac interoception (Attend, Feel, Regulate) while participants perform a ‘cardiac recognition’ task. 
For both the Feel and Regulate conditions, participants displayed enhanced recognition of their own heartbeat, 
accompanied by larger heartbeat-evoked potentials (HEPs), suggesting that these approaches could be used interchangeably. 
Importantly, meta-cognitive interoceptive insight was highest in the Regulate condition, indicative of stronger engagement with 
interoceptive signals in addition to greater ecological validity. Only in the passive interoception condition (Feel) was a 
significant association found between accuracy in recognising one’s own heartbeat and the amplitude of HEPs. 
Overall, our results imply that active conditions have an important role to play in future investigation of interoception.

Keywords: Interoception; Metacognition; Biofeedback; Predictive coding; Self-recognition

## Highlights

* We test cardioception when people attend, feel or regulate their HR.
* Feel and Regulate improve cardiac recognition.
* Regulate involves the highest levels of metacognition.
* Only in the Feel condition cardiac recognition correlates with HEPs.

## Coding of conditions 
* 1 = Attend, 2 = Feel, 3 = Regulate
* 1 = Congruent, 2 = Incongruent


## Description of Data and scripts
### HEP_diff_FvsA_dprime_diff_472484_AF8F4F6.xlsx
* To discover whether the heartbeat-evoked potential (HEP) amplitudes reflected behavioral differences, we investigated potential
links between cardiac recognition and the modulation of HEP amplitudes in each Condition. 
To match HEPs against d′– which inherently captures the Congruency to Incongruency relation –
we first calculated ‘Congruency Difference’ amplitude measures for HEPs in each of the three Conditions, by subtracting the mean 
amplitudes on Incongruent trials from those on  the  Congruent trials. Then, to  fully separate Condition-related effects from 
attentional processes, we  treated the Attend Condition as a baseline control (as it captured all the extero-ceptive aspects of the task)
and therefore subtracted the Congruency Difference amplitudes in the Attend condition from the Feel and Regulate Conditions (Fig. 5B).
To mirror this on a behavioral level, we subtracted d′scores in the Attend Condition from d′in the other two interoceptive Conditions respectively
(i.e., Feel and Regulate). 

### biosemi64_neighb.mat
* EEG layout map used in analysis.

### ecg_plot.m
* Script to create Figure 4B - Average ECG signal across all three Conditions (the solid line refers to the Congruent biofeedback and dashed lines 
to Incongruent feedback). Shaded areas around mean amplitudes indicate 95% confidence intervals. 

### indicesST4_wo_HEP_outliers.csv
* Signal detection indices per Engagement condition per participant.

### plot_eeg_dprime_for_regression.m
* The amplitude of Heartbeat Evoked Potential (HEP) difference measure was related to the difference in d′in the Feel compared with the Attend Condition.
Figure 5 A and B are plotting the relevant differences in amplitudes.

### plot_interaction.m
* Figure 4 A plotting the relevant differences in amplitudes under different Engagement strategies both for the Congruent and Incongruent conditions.

### plot_notincluded_groupedBYstrategies.m
* Like previous plot but signals are separated by Congruency condition then plotted for every Engagement condition (not included in final manuscript).

### study4_ECG_artifact_analysis.m
* Script used for ECG artifact analysis.

### study4_cbp_clean.m
* Script for the main EEG analysis.

### study4_eeg_v4.m
* The experiment's script.

### withall_CIC_all_avg_main_0405_Rfrontal.xlsx
* Amplitudes per condition for the right frontal ROI, within the a priori latency of 200–300 ms.

### Outliers
We used a multivariate model approach for outlier identification. Four influential outliers (ID = 5,11,16,30) were identified, based on the amount of
impact their data points had on the predicted outcome - represented by Cook’s distance (Cook, 1977). We decided to remove these participants as
they had more than one datapoint where Cook’s distance was four times greater than the mean, leaving us with a sample of N=30.


