%previous version eeg_pipeline_averaged.m
%in this analysis I use the ROI-s from Sel,Azevedo, Tsakiris(2017)
%Testing for main effects
%d-prime is regressed on synchrony difference score bw conditions same in
%Sel et al.(2017) F(C-IC)vs A(C-IC)
%Outliers were identified by inspecting the distribution in R using the
%raincloud plots - when Rfrontal are was inspected ID=5,11,16,30 were
%identified. They seem to be outliers in every conditions so difference
%scores should not be affected, yet from now on we exclude them from the
%analysis. For consistency these people
%will be excluded from all analyses.

clear all
ft_defaults %to avoid multiple versions
cd('D:\export2')

D = dir('*.dat');
numfiles=numel(D);
ALL_CA={};
ALL_CF={};
ALL_CR={};
ALL_ICA={};
ALL_ICF={};
ALL_ICR={};
%change this to exclude people, now excluded ID=5,11,16,30
sub=[1,2,3,4,6,7,8,9,10,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,31,32,33,34];

for  i=1:length(sub);
    q=sub(i)
    for j=1:6
cfg=[];
raw_filename=sprintf('%.0f_LH_Average_%d.dat',q,j); %just i i think if not using exclusion, for exclusion 'included(i)'
cfg.dataset=raw_filename;
cfg.channel = {'EEG', '-EXG*'}
eegData =ft_preprocessing(cfg);

% cfg.begsample = 50;
% cfg.endsample = 150;
% eegData_short = ft_redefinetrial(cfg,eegData);

ALL_S{j}= ft_timelockanalysis(cfg,eegData); %replace with eegData_short

    end
    ALL_CA{i}=ALL_S{1};
    ALL_CF{i}=ALL_S{2};
    ALL_CR{i}=ALL_S{3};
    ALL_ICA{i}=ALL_S{4};
    ALL_ICF{i}=ALL_S{5};
    ALL_ICR{i}=ALL_S{6};
end
 
ALL{1}=ALL_CA;
ALL{2}=ALL_CF;
ALL{3}=ALL_CR;
ALL{4}=ALL_ICA;
ALL{5}=ALL_ICF;
ALL{6}=ALL_ICR;

%to explore main effects
ALL_C=[ALL_CA ALL_CF ALL_CR];
ALL_IC=[ALL_ICA ALL_ICF ALL_ICR];
ALL_A= [ALL_CA ALL_ICA];
ALL_F= [ALL_CF ALL_ICF];
ALL_R= [ALL_CR ALL_ICR];

%to explore simple relationship between HEP and HRV and metacognition
ALL_COND_HEP=[ALL_CA ALL_CF ALL_CR ALL_ICA ALL_ICF ALL_ICR];
% save ALL_CA ALL_CA
% save ALL_CF ALL_CF
% save ALL_CR ALL_CR
% save ALL_ICA ALL_ICA
% save ALL_ICF ALL_ICF
% save ALL_ICR ALL_ICR
% save ALL_C ALL_C
% save ALL_IC ALL_IC
% save ALL_A ALL_A
% save ALL_F ALL_F
% save ALL_R ALL_R
% 
% load ALL_CA 
% load ALL_CF 
% load ALL_CR 
% load ALL_ICA 
% load ALL_ICF 
% load ALL_ICR 
% load ALL_C 
% load ALL_IC
% load ALL_A
% load ALL_F
% load ALL_R 

%6 rois sel et al
Lfrontal=	{'AF7', 'AF3', 'F5', 'F3'}
Lcentral= {'FC5', 'FC3', 'C5', 'C3'}
Lparietal=	{'CP5', 'CP3', 'P5', 'P3'}

Rfrontal= {'AF4', 'AF8', 'F4', 'F6'}
Rcentral= {'FC4', 'FC6', 'C4', 'C6'}
Rparietal= {'CP4', 'CP6', 'P4', 'P6'}

% Pollatos Schandry 2004 -> spread out activity but the activity on
% frontocentrasl areascould signal processing of cardio-afferent signals

% Lfrontocentral = {'FC1', 'FC5'}  %FC5: positive FC6: negative FC1:negative FC2:around 0 but positive
% Rfrontocentral = {'FC2', 'FC6'}
% frontocentral = {'FC1', 'FC5','FC2', 'FC6'}
load 'biosemi64_neighb.mat';
%%
% calculate the grand average for each condition
cfg = [];
%cfg.channel = frontal;
cfg.latency = [0 0.6]; %[0.3 0.6];%'all'; %CHANGE THIS WHEN CHANGING TIME WINDOW
cfg.parameter = 'avg';
GA_S1         = ft_timelockgrandaverage(cfg,ALL_CA{:});  
GA_S2        = ft_timelockgrandaverage(cfg,ALL_CF{:});
GA_S3        = ft_timelockgrandaverage(cfg,ALL_CR{:});
GA_S4        = ft_timelockgrandaverage(cfg,ALL_ICA{:});
GA_S5         = ft_timelockgrandaverage(cfg,ALL_ICF{:});  
GA_S6        = ft_timelockgrandaverage(cfg,ALL_ICR{:});

GA_C         = ft_timelockgrandaverage(cfg,ALL_C{:});
GA_IC        = ft_timelockgrandaverage(cfg,ALL_IC{:});
GA_A        = ft_timelockgrandaverage(cfg,ALL_A{:});
GA_F        = ft_timelockgrandaverage(cfg,ALL_F{:});
GA_R        = ft_timelockgrandaverage(cfg,ALL_R{:});

%%
% calculate the grand average for each condition for interaction
% with the field keepindividual set to ‘yes’.
cfg = [];
%cfg.channel = frontal;
%cfg.latency = [0 0.6]; %[0.3 0.6];%'all'; %CHANGE THIS WHEN CHANGING TIME WINDOW
%cfg.parameter = 'avg';
cfg.keepindividual = 'yes'
%cfg.method = 'within'
iGA_S1         = ft_timelockgrandaverage(cfg,ALL_CA{:});  
iGA_S2        = ft_timelockgrandaverage(cfg,ALL_CF{:});
iGA_S3        = ft_timelockgrandaverage(cfg,ALL_CR{:});
iGA_S4        = ft_timelockgrandaverage(cfg,ALL_ICA{:});
iGA_S5         = ft_timelockgrandaverage(cfg,ALL_ICF{:});  
iGA_S6        = ft_timelockgrandaverage(cfg,ALL_ICR{:});

iGA_A        = ft_timelockgrandaverage(cfg,ALL_A{:});
iGA_F        = ft_timelockgrandaverage(cfg,ALL_F{:});
iGA_R        = ft_timelockgrandaverage(cfg,ALL_R{:});
iGA_HEP      = ft_timelockgrandaverage(cfg,ALL_COND_HEP{:});

mean_ALL_COND = iGA_S1;

d_S2_S1=iGA_S2;
d_S3_S1=iGA_S3;
d_S5_S4=iGA_S5;
d_S6_S4=iGA_S6;
d_S3_S2=iGA_S3;
d_S6_S5=iGA_S6;

d_FvA=iGA_F;
d_RvA=iGA_F;

mean_ALL_COND.individual=(iGA_S1.individual+iGA_S2.individual+iGA_S3.individual+...
    iGA_S4.individual+iGA_S5.individual+iGA_S6.individual)/6;

d_S2_S1.individual=iGA_S2.individual-iGA_S1.individual;
d_S3_S1.individual=iGA_S3.individual-iGA_S1.individual;
d_S5_S4.individual=iGA_S5.individual-iGA_S4.individual;
d_S6_S4.individual=iGA_S6.individual-iGA_S4.individual;
d_S3_S2.individual=iGA_S3.individual-iGA_S2.individual;
d_S6_S5.individual=iGA_S6.individual-iGA_S5.individual;

d_FvA.individual=iGA_F.individual-iGA_A.individual;
d_RvA.individual=iGA_R.individual-iGA_A.individual;

d_S3S1_S2S1=iGA_S3;
d_S6S4_S5S4=iGA_S5;

d_S3S1_S2S1.individual=d_S3_S1.individual-d_S2_S1.individual;
d_S6S4_S5S4.individual=d_S6_S4.individual-d_S5_S4.individual;

% engagement effect
d_S1_S4=iGA_S1;
d_S2_S5=iGA_S2;
d_S3_S6=iGA_S3;

d_S1_S4.individual=iGA_S1.individual-iGA_S4.individual;
d_S2_S5.individual=iGA_S2.individual-iGA_S5.individual;
d_S3_S6.individual=iGA_S3.individual-iGA_S6.individual;

d_S2S5_S1S4=iGA_S2;
d_S3S6_S1S4=iGA_S3;

d_S2S5_S1S4.individual=d_S2_S5.individual-d_S1_S4.individual;
d_S3S6_S1S4.individual=d_S3_S6.individual-d_S1_S4.individual;

 
%% stat main effect CIC

cfg = [];
cfg.channel =  Rfrontal %frontal%%{'all'}
cfg.latency = [0.4 0.5]; %[0 1] or 'all' %0.3 to 0.5=100ms to 300ms after R
%cfg.method='analytic';
cfg.method = 'montecarlo';
cfg.avgovertime = 'no';
cfg.avgoverchan ='no';
cfg.statistic = 'ft_statfun_depsamplesT';  %'depsamplesFmultivariate'; for more than 2 conditions, for 2 depsamplesT, for analytic 'ft_statfun_depsamplesT'
%cfg.correctm = 'no';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = neighbours;  
cfg.tail = 0;%0 if comparing 2 conditions, 1 if more than 2
cfg.clustertail = 0; %0 if comparing 2 conditions, 1 if more than 2
cfg.alpha = 0.025;
% cfg.alpha = 0.05;
cfg.numrandomization = 10000;

subj = size(ALL_C,2); %or ALL_A
%for 2 conditions:
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

%%%
cfg.design = design; %for depsamplesF design3
cfg.uvar  = 1;
cfg.ivar  = 2;

%main effects
clear statCmain

[statCmain] = ft_timelockstatistics(cfg,ALL_C{:},ALL_IC{:}); %don't forget to change subject size

%% stat main effect AFR

cfg = [];
cfg.channel =  Rfrontal%{'AF4','F4','F6'} %frontal%%{'all'}
cfg.latency = [0.4 0.5]; %[0 1] or 'all' %0.3 to 0.5=100ms to 300ms after R
%cfg.method='analytic';
cfg.method = 'montecarlo';
cfg.avgovertime = 'no';
cfg.avgoverchan ='no';
cfg.statistic = 'ft_statfun_depsamplesT';  %'depsamplesFmultivariate'; for more than 2 conditions, for 2 depsamplesT, for analytic 'ft_statfun_depsamplesT'
%cfg.correctm = 'no';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = neighbours;  
cfg.tail = 0;%0 if comparing 2 conditions, 1 if more than 2
cfg.clustertail = 0; %0 if comparing 2 conditions, 1 if more than 2
cfg.alpha = 0.025;
% cfg.alpha = 0.05;
cfg.numrandomization = 10000;

subj = size(ALL_A,2); %or ALL_A
%for 2 conditions:
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

%%%
cfg.design = design; %for depsamplesF design3
cfg.uvar  = 1;
cfg.ivar  = 2;
clear statFAmain
clear statRAmain
clear statRFmain

[statFAmain] = ft_timelockstatistics(cfg,ALL_F{:},ALL_A{:}); %don't forget to change subject size
[statRAmain] = ft_timelockstatistics(cfg,ALL_R{:},ALL_A{:}); %don't forget to change subject size
[statRFmain] = ft_timelockstatistics(cfg,ALL_R{:},ALL_F{:}); %don't forget to change subject size

%diff from F and R as difference from A
clear statFRmain_diff
[statFRmain_diff] = ft_timelockstatistics(cfg,d_FvA, d_RvA); %don't forget to change subject size

%% as anova for AFR
%for 3 conditions
cfg = [];
cfg.channel =  Rfrontal %frontal%%{'all'}
cfg.latency = [0.4 0.5]; %[0 1] or 'all' %0.3 to 0.5=100ms to 300ms after R
%cfg.method='analytic';
cfg.method = 'montecarlo';
cfg.avgovertime = 'no';
cfg.avgoverchan ='no';
cfg.statistic = 'depsamplesFmultivariate';  %'depsamplesFmultivariate'; for more than 2 conditions, for 2 depsamplesT, for analytic 'ft_statfun_depsamplesT'
%cfg.correctm = 'no';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = neighbours;  
cfg.tail = 1;%0 if comparing 2 conditions, 1 if more than 2 %warning said:"For a dependent samples F-statistic, it does not make sense to calculate a two-sided critical value."
cfg.numrandomization = 10000;
cfg.clustertail = 1; %0 if comparing 2 conditions, 1 if more than 2
cfg.alpha = 0.05; 

subj = size(ALL_A,2); %or ALL_A
design = zeros(2,3*subj);
for i = 1:subj
  design(1,i) = i;
  design(1,subj+i) = i;
  design(1,(2*subj)+i) = i;
  %design(1,(3*subj)+i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
  design(1,(2*subj)+i) = i;
 % design(1,(3*subj)+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;
design(2,(2*subj)+1:3*subj) = 3;
%design(2,(3*subj)+1:4*subj) = 4;

cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

clear statAFRmain
[statAFRmain] = ft_timelockstatistics(cfg,ALL_A{:},ALL_F{:},ALL_R{:}); %significant

%% interaction stats
cfg = [];
cfg.channel =  Rfrontal % %%{'all'}
cfg.latency = [0.4 0.5]; %[0 1] or 'all' %0.3 to 0.5=100ms to 300ms after R
%cfg.method='analytic';
cfg.method = 'montecarlo';
cfg.avgovertime = 'no';
cfg.avgoverchan ='no';
cfg.statistic = 'ft_statfun_depsamplesT';  %'depsamplesFmultivariate'; for more than 2 conditions, for 2 depsamplesT, for analytic 'ft_statfun_depsamplesT'
%cfg.correctm = 'no';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = neighbours;  
cfg.tail = 0;%0 if comparing 2 conditions, 1 if more than 2
cfg.clustertail = 0; %0 if comparing 2 conditions, 1 if more than 2
cfg.alpha = 0.025;
% cfg.alpha = 0.05;
cfg.numrandomization = 10000;

subj = size(ALL_CF,2); %ALL_CA
%for 2 conditions:
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

%%%
cfg.design = design; %for depsamplesF design3
cfg.uvar  = 1;
cfg.ivar  = 2;
clear stat1
clear stat2
clear stat3
clear stat4a
clear stat5a
clear stat6a
clear stat7a
clear stat8a
clear stat9a

[stat1] = ft_timelockstatistics(cfg,  ALL_CA{:},ALL_ICA{:}); %
[stat2] = ft_timelockstatistics(cfg,  ALL_CF{:},ALL_ICF{:}); %
[stat3] = ft_timelockstatistics(cfg,  ALL_CR{:},ALL_ICR{:}); %
[stat4a] = ft_timelockstatistics(cfg,  ALL_CA{:},ALL_CF{:}); %
[stat5a] = ft_timelockstatistics(cfg,  ALL_CA{:},ALL_CR{:}); %
[stat6a] = ft_timelockstatistics(cfg,  ALL_CR{:},ALL_CF{:}); %
[stat7a] = ft_timelockstatistics(cfg,  ALL_ICA{:},ALL_ICF{:}); %
[stat8a] = ft_timelockstatistics(cfg,  ALL_ICA{:},ALL_ICR{:}); %
[stat9a] = ft_timelockstatistics(cfg,  ALL_ICR{:},ALL_ICF{:}); %

%% as anova for 6 conditions C/IC - AFR
%for 6 conditions
cfg = [];
cfg.channel =  Rfrontal %frontal%%{'all'}
cfg.latency = [0.4 0.5]; %[0 1] or 'all' %0.3 to 0.5=100ms to 300ms after R
%cfg.method='analytic';
cfg.method = 'montecarlo';
cfg.avgovertime = 'no';
cfg.avgoverchan ='no';
cfg.statistic = 'depsamplesFmultivariate';  %'depsamplesFmultivariate'; for more than 2 conditions, for 2 depsamplesT, for analytic 'ft_statfun_depsamplesT'
%cfg.correctm = 'no';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.neighbours = neighbours;  
cfg.tail = 1;%0 if comparing 2 conditions, 1 if more than 2 %warning said:"For a dependent samples F-statistic, it does not make sense to calculate a two-sided critical value."
cfg.numrandomization = 20000;
cfg.clustertail = 1; %0 if comparing 2 conditions, 1 if more than 2
cfg.alpha = 0.05; 

subj = size(ALL_CA,2); %or ALL_A
design = zeros(2,3*subj);
for i = 1:subj
  design(1,i) = i;
  design(1,subj+i) = i;
  design(1,(2*subj)+i) = i;
  design(1,(3*subj)+i) = i;
  design(1,(4*subj)+i) = i;
  design(1,(5*subj)+i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
  design(1,(2*subj)+i) = i;
  design(1,(3*subj)+i) = i;
  design(1,(4*subj)+i) = i;
  design(1,(5*subj)+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;
design(2,(2*subj)+1:3*subj) = 3;
design(2,(3*subj)+1:4*subj) = 4;
design(2,(4*subj)+1:5*subj) = 5;
design(2,(5*subj)+1:6*subj) = 6;
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

clear statCICAFRmain
[statCICAFRmain] = ft_timelockstatistics(cfg,ALL_CA{:},ALL_CF{:},ALL_CR{:},ALL_ICA{:},ALL_ICF{:},ALL_ICR{:});

%% regression with d prime
% prepare continious ivar-s
indicesST4 = importdata('indicesST4_wo_HEP_outliers.csv'); %indicesST4_wo_HEP_outliers.csv
att=find(indicesST4.data(:,end)==1);
feel=find(indicesST4.data(:,end)==2);
reg=find(indicesST4.data(:,end)==3);

dprime_att=[indicesST4.data(att)];
dprime_feel=[indicesST4.data(feel)];
dprime_reg=[indicesST4.data(reg)];

cfg = [];
cfg.channel = Rfrontal
cfg.latency = [0.432 0.480]; %[0 1] or 'all' %0.3 to 0.5=100ms to 300ms after R
cfg.method = 'montecarlo';
cfg.avgovertime='no';
cfg.avgoverchan ='no';
cfg.statistic = 'ft_statfun_indepsamplesregrT'; %'depsamplesFmultivariate'; for more than 2 conditions, for 2 depsamplesT, regression 'indepsamplesregrT'
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 1;
cfg.neighbours = neighbours;  
cfg.tail = 0;%0 if comparing 2 conditions, 1 if more than 2
cfg.clustertail = 0; %0 if comparing 2 conditions, 1 if more than 2
cfg.alpha = 0.025; %positive relationship expected
cfg.numrandomization = 15000;
subj = size(ALL_CF,2);

design = zeros(2,subj);
for i = 1:subj
  design(1,i) = i;
end
cfg.ivar  = 2;

clear statD1
clear statD2
clear statD3
clear statD4
clear statD5

design(2,1:subj)   = dprime_att;
cfg.design = design;
[statD1] = ft_timelockstatistics(cfg,  d_S1_S4); %use dprime_att above

design(2,1:subj)   = dprime_feel;
cfg.design = design;
[statD2] = ft_timelockstatistics(cfg,  d_S2_S5); %use dprime_feel above

design(2,1:subj)   = dprime_reg;
cfg.design = design;
[statD3] = ft_timelockstatistics(cfg,  d_S3_S6); %use dprime_reg above

design(2,1:subj)   = dprime_feel-dprime_att;
cfg.design = design;
[statD4] = ft_timelockstatistics(cfg,  d_S2S5_S1S4); %equivalent of Sel et al.(2017)

design(2,1:subj)   = dprime_reg-dprime_att;
cfg.design = design;
[statD5] = ft_timelockstatistics(cfg,  d_S3S6_S1S4); %equivalent of Sel et al.(2017)


%% to get clusters and their time window
%just change stat to the stat interest
stat_name={'statCmain','statFAmain','statRAmain','statRFmain',...
            'stat1','stat2','stat3','stat4a','stat5a','stat6a','stat7a','stat8a','stat9a',...
            'statD1','statD2','statD3','statD4','statD5'};
stat_all={statCmain,statFAmain,statRAmain,statRFmain,...
            stat1,stat2,stat3,stat4a,stat5a,stat6a,stat7a,stat8a,stat9a,...
            statD1,statD2,statD3,statD4,statD5};        
stat_summary={};
for  i=1:length(stat_all);
pos_cluster_pvals={};
pos_signif_clust={};
pos_Tsum={};
pos={};
neg_cluster_pvals={};
neg_signif_clust={};
neg_Tsum={};
neg={};
clust_chan_pos={};
clust_chan_neg={};
clust_timewindow_pos={};
clust_timewindow_neg={};

stat= stat_all{i};
if length(fieldnames(stat))>= 15
%positive cluster stats
if isempty(stat.posclusters(:))==0
pos_cluster_pvals{1} = [stat.posclusters(:).prob];
pos_Tsum{1} = [stat.posclusters(:).clusterstat];
pos_signif_clust{1} = find(pos_cluster_pvals{1} < stat.cfg.alpha);
pos{1} = ismember(stat.posclusterslabelmat, pos_signif_clust{1});
[clust_chan,ms]=find(stat.posclusterslabelmat==1);
clust_chan=unique(clust_chan);
clust_timewindow_pos=[min(unique(ms)) max(unique(ms))];
clust_chan_pos{1}=stat.label(clust_chan).';
else
    pos_cluster_pvals{1}='Nan';
    pos_Tsum{1} = 'Nan';
    pos_signif_clust{1}='Nan';
    clust_chan_pos{1}='Nan';
    clust_timewindow_pos{1}='Nan';
end

%negative cluster stats
if isempty(stat.negclusters(:))==0
neg_cluster_pvals{1} = [stat.negclusters(:).prob];
neg_Tsum{1} = [stat.negclusters(:).clusterstat];
neg_signif_clust{1} = find(neg_cluster_pvals{1} < stat.cfg.alpha);
neg{1} = ismember(stat.negclusterslabelmat, neg_signif_clust{1});

[clust_chan,ms]=find(stat.negclusterslabelmat==1);
clust_chan=unique(clust_chan);
clust_timewindow_neg=[min(unique(ms)) max(unique(ms))];
clust_chan_neg{1}=stat.label(clust_chan).';
else
    neg_cluster_pvals{1}='Nan';
    neg_Tsum{1} = 'Nan';
    neg_signif_clust{1}='Nan';
    clust_chan_neg{1}='Nan';
    clust_timewindow_neg{1}='Nan';
end
    
stat_summary{i,1}=stat_name{i};
stat_summary{i,2}=pos_cluster_pvals{:};
stat_summary{i,3}=pos_Tsum{:}
stat_summary{i,4}=pos_signif_clust{:};
stat_summary{i,5}=clust_chan_pos{:};
stat_summary{i,6}=clust_timewindow_pos(:);

stat_summary{i,7}=neg_cluster_pvals{:};
stat_summary{i,8}=neg_Tsum{:};
stat_summary{i,9}=neg_signif_clust{:};
stat_summary{i,10}=clust_chan_neg{:};
stat_summary{i,11}=clust_timewindow_neg(:);
  
else
stat_summary{i,1}=stat_name{i};
stat_summary{i,2}='Nan';
stat_summary{i,3}='Nan';

stat_summary{i,4}='Nan';
stat_summary{i,5}='Nan';
stat_summary{i,6}='Nan';

stat_summary{i,7}='Nan';
stat_summary{i,8}='Nan';

stat_summary{i,9}='Nan';
stat_summary{i,10}='Nan';
stat_summary{i,11}='Nan';
end
end

stat_name2={'statAFRmain','statCICAFRmain'};
stat_all2={statAFRmain,statCICAFRmain};        
for  v=1:length(stat_all2);
pos_cluster_pvals={};
pos_Tsum={};
pos_signif_clust={};
pos={};
clust_chan_pos={};
clust_timewindow_pos={};

stat= stat_all2{v};
pos_cluster_pvals{1} = [stat.posclusters(:).prob];
pos_Tsum{1} = [stat.posclusters(:).clusterstat];
pos_signif_clust{1} = find(pos_cluster_pvals{1} < stat.cfg.alpha);
pos{1} = ismember(stat.posclusterslabelmat, pos_signif_clust{1});
[clust_chan,ms]=find(stat.posclusterslabelmat==1);
clust_chan=unique(clust_chan);
clust_timewindow_pos=[min(unique(ms)) max(unique(ms))];
clust_chan_pos{1}=stat.label(clust_chan).';

stat_summary{i+v,1}=stat_name2{v};
stat_summary{i+v,2}=pos_cluster_pvals{:};
stat_summary{i+v,3}=pos_Tsum{:};
stat_summary{i+v,4}=pos_signif_clust{:};
stat_summary{i+v,5}=clust_chan_pos{:};
stat_summary{i+v,6}=clust_timewindow_pos(:);

end
save stat_summary stat_summary
writetable(table(stat_summary),'stat_fieldtrip_summary_study4.xlsx','WriteVariableNames',true) %with varnames as header


%% create spreadsheet for R
%all 2x3 conditions
for ptp=1:size(ALL_C,2)
    
cfg = [];
%cfg.channel = n22_new_fearanger_common;%clust_chan;
cfg.latency = [0.4 0.5]; %[0.3 0.6];%'all'; %CHANGE THIS WHEN CHANGING TIME WINDOW
cfg.parameter = 'avg';
cfg.channel = Rfrontal;

avg_C(ptp)       = ft_timelockgrandaverage(cfg,ALL_C{ptp});
avg_IC(ptp)       = ft_timelockgrandaverage(cfg,ALL_IC{ptp});

ptp_avg_C(ptp)=mean(mean(avg_C(ptp).avg,1));
ptp_avg_IC(ptp)=mean(mean(avg_IC(ptp).avg,1));
 
end

ptp_avg_C=ptp_avg_C.';
ptp_avg_IC=ptp_avg_IC.';

ID=transpose([1:34,1:34,1:34,1:34,1:34,1:34]);
%1 c 2% ic
Congruency=[ones(3*34,1);2*ones(3*34,1)];
%1 a 2 f 3 r
Strategy= [ones(34,1);2*ones(34,1);3*ones(34,1);1*ones(34,1);2*ones(34,1);3*ones(34,1)];

CIC_all_avg_main=table([ptp_avg_C;ptp_avg_IC],ID,Congruency,Strategy);

filename = 'CIC_all_avg_main_0405_Rfrontal.xlsx';
writetable(CIC_all_avg_main,filename) 
%%
%calculating averages for R
clear avg_S1
clear avg_S2
clear avg_S3
clear avg_S4
clear avg_S5
clear avg_S6

clear ptp_avg_S1
clear ptp_avg_S2
clear ptp_avg_S3
clear ptp_avg_S4
clear ptp_avg_S5
clear ptp_avg_S6

%d prime
for ptp=1:size(ALL_CA,2)
    
cfg = [];
cfg.latency =[0.472 0.484]%
cfg.parameter = 'avg';
cfg.channel =  {'AF8' 'F4' 'F6'}
avg_S1(ptp)   = ft_timelockgrandaverage(cfg,ALL_CA{ptp});%{baby}  
avg_S2(ptp)      = ft_timelockgrandaverage(cfg,ALL_CF{ptp});
avg_S3(ptp)        = ft_timelockgrandaverage(cfg,ALL_CR{ptp});
avg_S4(ptp)        = ft_timelockgrandaverage(cfg,ALL_ICA{ptp});
avg_S5(ptp)         = ft_timelockgrandaverage(cfg,ALL_ICF{ptp});  
avg_S6(ptp)       = ft_timelockgrandaverage(cfg,ALL_ICR{ptp});


    
 ptp_avg_S1(ptp)=mean(mean(avg_S1(ptp).avg,1));
 ptp_avg_S2(ptp)=mean(mean(avg_S2(ptp).avg,1));
 ptp_avg_S3(ptp)=mean(mean(avg_S3(ptp).avg,1));
 ptp_avg_S4(ptp)=mean(mean(avg_S4(ptp).avg,1));
 ptp_avg_S5(ptp)=mean(mean(avg_S5(ptp).avg,1));
 ptp_avg_S6(ptp)=mean(mean(avg_S6(ptp).avg,1));
 
end
ptp_avg_S1=ptp_avg_S1.';
ptp_avg_S2=ptp_avg_S2.';
ptp_avg_S3=ptp_avg_S3.';
ptp_avg_S4=ptp_avg_S4.';
ptp_avg_S5=ptp_avg_S5.';
ptp_avg_S6=ptp_avg_S6.';

ptp_avg_S2S5vsS1S4 = (ptp_avg_S2-ptp_avg_S5)-(ptp_avg_S1-ptp_avg_S4);
ptp_avg_S3S6vsS1S4 = (ptp_avg_S3-ptp_avg_S6)-(ptp_avg_S1-ptp_avg_S4);

diff_dprime_FA=dprime_feel-dprime_att;
ID=transpose(sub);
Dprime_all_avg_main=table(ID,ptp_avg_S2S5vsS1S4,diff_dprime_FA);
filename = 'HEP_diff_FvsA_dprime_diff_472484_AF8F4F6.xlsx';
writetable(Dprime_all_avg_main,filename)

%% figures
%plot singlpe plot for erp-s and topoplot

figure
hold on
cfg = [];
cfg.xlim = [0.1 1];
%cfg.ylim = [-0.025 -0.005];
cfg.graphcolor    ='bb';
cfg.linewidth = 1.5;
cfg.linestyle     = {'-' '--'};
cfg.layout='biosemi64.lay';
cfg.channel = Rfrontal;
colorbar
clf;
ft_singleplotER(cfg,GA_C,GA_IC); %,GA_S4,GA_S5,GA_S6

figure
hold on
cfg = [];
cfg.xlim = [0.1 1];
%cfg.ylim = [-0.025 -0.005];
cfg.graphcolor    ='gbr';
cfg.linewidth = 1.5;
cfg.linestyle     = {'-' '-' '-'};
cfg.layout='biosemi64.lay';
cfg.channel = Rfrontal;
colorbar
clf;
ft_singleplotER(cfg,GA_A,GA_F,GA_R); %,GA_S4,GA_S5,GA_S6


figure
hold on
cfg = [];
cfg.xlim = [0.1 1];
cfg.ylim = [-2 1];
cfg.graphcolor    ='gbrgbr';
cfg.linewidth = 1.5;
cfg.linestyle     = {'-' '-' '-' '--' '--' '--'};
cfg.layout='biosemi64.lay';
cfg.channel = Rfrontal;
colorbar
clf;
ft_singleplotER(cfg,GA_S1,GA_S2,GA_S3,GA_S4,GA_S5,GA_S6); %,GA_S4,GA_

%prepare for plotting
%first level
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
%prep
GA_S2vsS1 = ft_math(cfg,GA_S2, GA_S1);
GA_S3vsS1 = ft_math(cfg,GA_S3, GA_S1);
GA_S3vsS2 = ft_math(cfg,GA_S3, GA_S2);

GA_S5vsS4 = ft_math(cfg,GA_S5, GA_S4);
GA_S6vsS4 = ft_math(cfg,GA_S6, GA_S4);
GA_S6vsS5 = ft_math(cfg,GA_S6, GA_S5);

%%same but flipped
GA_S4vsS5 = ft_math(cfg,GA_S4, GA_S5);
GA_S4vsS6 = ft_math(cfg,GA_S4, GA_S6);

GA_S1vsS4 = ft_math(cfg,GA_S1, GA_S4);
GA_S3vsS6 = ft_math(cfg,GA_S3, GA_S6);
GA_S2vsS5 = ft_math(cfg,GA_S2, GA_S5);

GA_FvsA = ft_math(cfg,GA_F, GA_A);
GA_RvsA = ft_math(cfg,GA_R, GA_A);
GA_FvsR = ft_math(cfg,GA_F, GA_R);


%second level
%stat1
GA_S2S1vsS5S4 = ft_math(cfg,GA_S2vsS1,GA_S5vsS4);
%stat2
GA_S3S1vsS6S4 = ft_math(cfg,GA_S3vsS1,GA_S6vsS4);
%stat3
GA_S3S2vsS6S5 = ft_math(cfg,GA_S3vsS2,GA_S6vsS5);
%stat4
GA_S1S4vsS2S5 = ft_math(cfg,GA_S1vsS4,GA_S2vsS5);
%stat5
GA_S1S4vsS3S6 = ft_math(cfg,GA_S1vsS4,GA_S3vsS6);
%stat6
GA_S3S6vsS2S5 = ft_math(cfg,GA_S3vsS6,GA_S2vsS5);
%prep
GA_S2S5vsS1S4 = ft_math(cfg,GA_S2vsS5,GA_S1vsS4);
GA_S3S6vsS1S4 = ft_math(cfg,GA_S3vsS6,GA_S1vsS4);

%third level
%stat7
GA_S2S5S1S4vsS3S6S1S4 = ft_math(cfg,GA_S2S5vsS1S4,GA_S3S6vsS1S4);


%from tutorial %works %chnage parameters according to interest
cfg = [];
cfg.xlim = [0.4 0.5]; %0.4 0.5;0.4 0.468;0.4 0.452
cfg.zlim = [-0.5 0.5];
%cfg.channel = Rparietal;
cfg.layout = 'biosemi64.lay';
cfg.parameter = 'avg'; % the default 'avg' 'individual' is not present in the data
figure; ft_topoplotER(cfg,GA_S3S6vsS1S4); colorbar%GA_S2S5vsS1S4 %GA_S3S6vsS1S4

% so far it was the same as above, now change the colormap
ft_hastoolbox('brewermap', 1); % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

figure
hold on
cfg = [];
cfg.xlim = [0 0.6];
cfg.ylim = [-0.5 0.5];
cfg.graphcolor    ='gbr';
cfg.linewidth = 1.5;
cfg.linestyle     = {'-' '-' '-'};
cfg.layout='biosemi64.lay';
cfg.channel = Rfrontal;
colorbar
clf;
ft_singleplotER(cfg,GA_S1vsS4,GA_S2vsS5,GA_S3vsS6); %,GA_S4,GA_
