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
ALL_CA_ECG={};
ALL_CF_ECG={};
ALL_CR_ECG={};
ALL_ICA_ECG={};
ALL_ICF_ECG={};
ALL_ICR_ECG={};
%change this to exclude people, now excluded ID=5,11,16,30
sub=[1,2,3,4,6,7,8,9,10,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,29,31,32,33,34];

for  i=1:length(sub);
    q=sub(i)
    for j=1:6
cfg=[];
raw_filename=sprintf('%.0f_LH_Average_%d.dat',q,j); %just i i think if not using exclusion, for exclusion 'included(i)'
cfg.dataset=raw_filename;
cfg.channel = {'-EEG', 'EXG8'}
eegData =ft_preprocessing(cfg);



ALL_S_ECG{j}= ft_timelockanalysis(cfg,eegData); %replace with eegData_short

    end
    ALL_CA_ECG{i}=ALL_S_ECG{1};
    ALL_CF_ECG{i}=ALL_S_ECG{2};
    ALL_CR_ECG{i}=ALL_S_ECG{3};
    ALL_ICA_ECG{i}=ALL_S_ECG{4};
    ALL_ICF_ECG{i}=ALL_S_ECG{5};
    ALL_ICR_ECG{i}=ALL_S_ECG{6};
end
 
ALL{1}=ALL_CA_ECG;
ALL{2}=ALL_CF_ECG;
ALL{3}=ALL_CR_ECG;
ALL{4}=ALL_ICA_ECG;
ALL{5}=ALL_ICF_ECG;
ALL{6}=ALL_ICR_ECG;

%to explore main effects
ALL_C_ECG=[ALL_CA_ECG ALL_CF_ECG ALL_CR_ECG];
ALL_IC_ECG=[ALL_ICA_ECG ALL_ICF_ECG ALL_ICR_ECG];
ALL_A_ECG= [ALL_CA_ECG ALL_ICA_ECG];
ALL_F_ECG= [ALL_CF_ECG ALL_ICF_ECG];
ALL_R_ECG= [ALL_CR_ECG ALL_ICR_ECG];

%to explore simple relationship between HEP and HRV and metacognition
ALL_COND_HEP_ECG=[ALL_CA_ECG ALL_CF_ECG ALL_CR_ECG ALL_ICA_ECG ALL_ICF_ECG ALL_ICR_ECG];

%roi for cardiac artefact
ECG = {'EXG8'}
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
GA_S1_ECG         = ft_timelockgrandaverage(cfg,ALL_CA_ECG{:});  
GA_S2_ECG        = ft_timelockgrandaverage(cfg,ALL_CF_ECG{:});
GA_S3_ECG        = ft_timelockgrandaverage(cfg,ALL_CR_ECG{:});
GA_S4_ECG        = ft_timelockgrandaverage(cfg,ALL_ICA_ECG{:});
GA_S5_ECG         = ft_timelockgrandaverage(cfg,ALL_ICF_ECG{:});  
GA_S6_ECG        = ft_timelockgrandaverage(cfg,ALL_ICR_ECG{:});

GA_C_ECG         = ft_timelockgrandaverage(cfg,ALL_C_ECG{:});
GA_IC_ECG        = ft_timelockgrandaverage(cfg,ALL_IC_ECG{:});
GA_A_ECG        = ft_timelockgrandaverage(cfg,ALL_A_ECG{:});
GA_F_ECG        = ft_timelockgrandaverage(cfg,ALL_F_ECG{:});
GA_R_ECG        = ft_timelockgrandaverage(cfg,ALL_R_ECG{:});

%%
% calculate the grand average for each condition for interaction
% with the field keepindividual set to ‘yes’.
cfg = [];
cfg.keepindividual = 'yes'
iGA_S1_ECG         = ft_timelockgrandaverage(cfg,ALL_CA_ECG{:});  
iGA_S2_ECG        = ft_timelockgrandaverage(cfg,ALL_CF_ECG{:});
iGA_S3_ECG        = ft_timelockgrandaverage(cfg,ALL_CR_ECG{:});
iGA_S4_ECG        = ft_timelockgrandaverage(cfg,ALL_ICA_ECG{:});
iGA_S5_ECG         = ft_timelockgrandaverage(cfg,ALL_ICF_ECG{:});  
iGA_S6_ECG        = ft_timelockgrandaverage(cfg,ALL_ICR_ECG{:});

iGA_A_ECG        = ft_timelockgrandaverage(cfg,ALL_A_ECG{:});
iGA_F_ECG        = ft_timelockgrandaverage(cfg,ALL_F_ECG{:});
iGA_R_ECG        = ft_timelockgrandaverage(cfg,ALL_R_ECG{:});
iGA_HEP_ECG      = ft_timelockgrandaverage(cfg,ALL_COND_HEP_ECG{:});

mean_ALL_COND_ECG = iGA_S1_ECG;

d_S2_S1_ECG=iGA_S2_ECG;
d_S3_S1_ECG=iGA_S3_ECG;
d_S5_S4_ECG=iGA_S5_ECG;
d_S6_S4_ECG=iGA_S6_ECG;
d_S3_S2_ECG=iGA_S3_ECG;
d_S6_S5_ECG=iGA_S6_ECG;

d_FvA_ECG=iGA_F_ECG;
d_RvA_ECG=iGA_F_ECG;

mean_ALL_COND_ECG.individual=(iGA_S1_ECG.individual+iGA_S2_ECG.individual+iGA_S3_ECG.individual+...
    iGA_S4_ECG.individual+iGA_S5_ECG.individual+iGA_S6_ECG.individual)/6;

d_S2_S1_ECG.individual=iGA_S2_ECG.individual-iGA_S1_ECG.individual;
d_S3_S1_ECG.individual=iGA_S3_ECG.individual-iGA_S1_ECG.individual;
d_S5_S4_ECG.individual=iGA_S5_ECG.individual-iGA_S4_ECG.individual;
d_S6_S4_ECG.individual=iGA_S6_ECG.individual-iGA_S4_ECG.individual;
d_S3_S2_ECG.individual=iGA_S3_ECG.individual-iGA_S2_ECG.individual;
d_S6_S5_ECG.individual=iGA_S6_ECG.individual-iGA_S5_ECG.individual;

d_FvA_ECG.individual=iGA_F_ECG.individual-iGA_A_ECG.individual;
d_RvA_ECG.individual=iGA_R_ECG.individual-iGA_A_ECG.individual;

d_S3S1_S2S1_ECG=iGA_S3_ECG;
d_S6S4_S5S4_ECG=iGA_S5_ECG;

d_S3S1_S2S1_ECG.individual=d_S3_S1_ECG.individual-d_S2_S1_ECG.individual;
d_S6S4_S5S4_ECG.individual=d_S6_S4_ECG.individual-d_S5_S4_ECG.individual;

% engagement effect
d_S1_S4_ECG=iGA_S1_ECG;
d_S2_S5_ECG=iGA_S2_ECG;
d_S3_S6_ECG=iGA_S3_ECG;

d_S1_S4_ECG.individual=iGA_S1_ECG.individual-iGA_S4_ECG.individual;
d_S2_S5_ECG.individual=iGA_S2_ECG.individual-iGA_S5_ECG.individual;
d_S3_S6_ECG.individual=iGA_S3_ECG.individual-iGA_S6_ECG.individual;

d_S2S5_S1S4_ECG=iGA_S2_ECG;
d_S3S6_S1S4_ECG=iGA_S3_ECG;

d_S2S5_S1S4_ECG.individual=d_S2_S5_ECG.individual-d_S1_S4_ECG.individual;
d_S3S6_S1S4_ECG.individual=d_S3_S6_ECG.individual-d_S1_S4_ECG.individual;


%% as anova for 6 conditions C/IC - AFR
%for 6 conditions
cfg = [];
cfg.channel =  ECG 
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

subj = size(ALL_CA_ECG,2); %or ALL_A
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

clear statCICAFRmain_ECG
[statCICAFRmain_ECG] = ft_timelockstatistics(cfg,ALL_CA_ECG{:},ALL_CF_ECG{:},ALL_CR_ECG{:},ALL_ICA_ECG{:},ALL_ICF_ECG{:},ALL_ICR_ECG{:});%no effect

% %% regression with d prime
% % prepare continious ivar-s
% indicesST4 = importdata('indicesST4_wo_HEP_outliers.csv'); %indicesST4_wo_HEP_outliers.csv
% att=find(indicesST4.data(:,end)==1);
% feel=find(indicesST4.data(:,end)==2);
% reg=find(indicesST4.data(:,end)==3);
% 
% dprime_att=[indicesST4.data(att)];
% dprime_feel=[indicesST4.data(feel)];
% dprime_reg=[indicesST4.data(reg)];
% 
% %difference reg
% 
% cfg = [];
% cfg.channel = ECG
% cfg.latency = [0.432 0.480]; %[0 1] or 'all' %0.3 to 0.5=100ms to 300ms after R
% cfg.method = 'montecarlo';
% cfg.avgovertime='no';
% cfg.avgoverchan ='no';
% cfg.statistic = 'ft_statfun_indepsamplesregrT'; %'depsamplesFmultivariate'; for more than 2 conditions, for 2 depsamplesT, regression 'indepsamplesregrT'
% cfg.correctm = 'cluster';
% cfg.clusteralpha = 0.05;
% cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan = 1;
% cfg.neighbours = neighbours;  
% cfg.tail = 0;%0 if comparing 2 conditions, 1 if more than 2
% cfg.clustertail = 0; %0 if comparing 2 conditions, 1 if more than 2
% cfg.alpha = 0.025; %positive relationship expected
% cfg.numrandomization = 15000;
% subj = size(ALL_CF_ECG,2); %d_S2_S5 %d_S1_S4
% 
% design = zeros(2,subj);
% for i = 1:subj
%   design(1,i) = i;
% end
% cfg.ivar  = 2;
% 
% clear statD4_ECG
% 
% design(2,1:subj)   = dprime_feel-dprime_att;
% cfg.design = design;
% [statD4_ECG] = ft_timelockstatistics(cfg,  d_S2S5_S1S4_ECG); %equivalent of Sel et al.(2017)


%%
% %calculating averages for R
% clear avg_S1_ECG
% clear avg_S2_ECG
% clear avg_S3_ECG
% clear avg_S4_ECG
% clear avg_S5_ECG
% clear avg_S6_ECG
% 
% clear ptp_avg_S1_ECG
% clear ptp_avg_S2_ECG
% clear ptp_avg_S3_ECG
% clear ptp_avg_S4_ECG
% clear ptp_avg_S5_ECG
% clear ptp_avg_S6_ECG
% 
% for ptp=1:size(ALL_CA_ECG,2)
%     
% cfg = [];
% cfg.latency =[0.472 0.484]%
% cfg.parameter = 'avg';
% cfg.channel =  ECG
% avg_S1_ECG(ptp)   = ft_timelockgrandaverage(cfg,ALL_CA_ECG{ptp});%{baby}  
% avg_S2_ECG(ptp)      = ft_timelockgrandaverage(cfg,ALL_CF_ECG{ptp});
% avg_S3_ECG(ptp)        = ft_timelockgrandaverage(cfg,ALL_CR_ECG{ptp});
% avg_S4_ECG(ptp)        = ft_timelockgrandaverage(cfg,ALL_ICA_ECG{ptp});
% avg_S5_ECG(ptp)         = ft_timelockgrandaverage(cfg,ALL_ICF_ECG{ptp});  
% avg_S6_ECG(ptp)       = ft_timelockgrandaverage(cfg,ALL_ICR_ECG{ptp});
% 
% 
%     
%  ptp_avg_S1_ECG(ptp)=mean(mean(avg_S1_ECG(ptp).avg,1));
%  ptp_avg_S2_ECG(ptp)=mean(mean(avg_S2_ECG(ptp).avg,1));
%  ptp_avg_S3_ECG(ptp)=mean(mean(avg_S3_ECG(ptp).avg,1));
%  ptp_avg_S4_ECG(ptp)=mean(mean(avg_S4_ECG(ptp).avg,1));
%  ptp_avg_S5_ECG(ptp)=mean(mean(avg_S5_ECG(ptp).avg,1));
%  ptp_avg_S6_ECG(ptp)=mean(mean(avg_S6_ECG(ptp).avg,1));
%  
% end
% ptp_avg_S1_ECG=ptp_avg_S1_ECG.';
% ptp_avg_S2_ECG=ptp_avg_S2_ECG.';
% ptp_avg_S3_ECG=ptp_avg_S3_ECG.';
% ptp_avg_S4_ECG=ptp_avg_S4_ECG.';
% ptp_avg_S5_ECG=ptp_avg_S5_ECG.';
% ptp_avg_S6_ECG=ptp_avg_S6_ECG.';
% 
% ptp_avg_S2S5vsS1S4_ECG = (ptp_avg_S2_ECG-ptp_avg_S5_ECG)-(ptp_avg_S1_ECG-ptp_avg_S4_ECG);
% ptp_avg_S3S6vsS1S4_ECG = (ptp_avg_S3_ECG-ptp_avg_S6_ECG)-(ptp_avg_S1_ECG-ptp_avg_S4_ECG);
% 
% diff_dprime_FA_ECG=dprime_feel-dprime_att;
% ID=transpose(sub);
% Dprime_all_avg_main_ECG=table(ID,ptp_avg_S2S5vsS1S4_ECG,diff_dprime_FA_ECG);
% filename = 'ECG_diff_FvsA_dprime_diff_472484_AF8F4F6.xlsx';
% writetable(Dprime_all_avg_main_ECG,filename)


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
ft_singleplotER(cfg,GA_C_ECG,GA_IC_ECG); %,GA_S4,GA_S5,GA_S6

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
ft_singleplotER(cfg,GA_A_ECG,GA_F_ECG,GA_R_ECG); %,GA_S4,GA_S5,GA_S6


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
ft_singleplotER(cfg,GA_S1_ECG,GA_S2_ECG,GA_S3_ECG,GA_S4_ECG,GA_S5_ECG,GA_S6_ECG); %,GA_S4,GA_

%prepare for plotting
%first level
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
%prep
GA_S2vsS1_ECG = ft_math(cfg,GA_S2_ECG, GA_S1_ECG);
GA_S3vsS1_ECG = ft_math(cfg,GA_S3_ECG, GA_S1_ECG);
GA_S3vsS2_ECG = ft_math(cfg,GA_S3_ECG, GA_S2_ECG);

GA_S5vsS4_ECG = ft_math(cfg,GA_S5_ECG, GA_S4_ECG);
GA_S6vsS4_ECG = ft_math(cfg,GA_S6_ECG, GA_S4_ECG);
GA_S6vsS5_ECG = ft_math(cfg,GA_S6_ECG, GA_S5_ECG);

%%same but flipped
GA_S4vsS5_ECG = ft_math(cfg,GA_S4_ECG, GA_S5_ECG);
GA_S4vsS6_ECG = ft_math(cfg,GA_S4_ECG, GA_S6_ECG);

GA_S1vsS4_ECG = ft_math(cfg,GA_S1_ECG, GA_S4_ECG);
GA_S3vsS6_ECG = ft_math(cfg,GA_S3_ECG, GA_S6_ECG);
GA_S2vsS5_ECG = ft_math(cfg,GA_S2_ECG, GA_S5_ECG);

GA_FvsA_ECG = ft_math(cfg,GA_F_ECG, GA_A_ECG);
GA_RvsA_ECG = ft_math(cfg,GA_R_ECG, GA_A_ECG);
GA_FvsR_ECG = ft_math(cfg,GA_F_ECG, GA_R_ECG);


%second level
%stat1
GA_S2S1vsS5S4_ECG = ft_math(cfg,GA_S2vsS1_ECG,GA_S5vsS4_ECG);
%stat2
GA_S3S1vsS6S4_ECG = ft_math(cfg,GA_S3vsS1_ECG,GA_S6vsS4_ECG);
%stat3
GA_S3S2vsS6S5_ECG = ft_math(cfg,GA_S3vsS2_ECG,GA_S6vsS5_ECG);
%stat4
GA_S1S4vsS2S5_ECG = ft_math(cfg,GA_S1vsS4_ECG,GA_S2vsS5_ECG);
%stat5
GA_S1S4vsS3S6_ECG = ft_math(cfg,GA_S1vsS4_ECG,GA_S3vsS6_ECG);
%stat6
GA_S3S6vsS2S5_ECG = ft_math(cfg,GA_S3vsS6_ECG,GA_S2vsS5_ECG);
%prep
GA_S2S5vsS1S4_ECG = ft_math(cfg,GA_S2vsS5_ECG,GA_S1vsS4_ECG);
GA_S3S6vsS1S4_ECG = ft_math(cfg,GA_S3vsS6_ECG,GA_S1vsS4_ECG);

%third level
%stat7
GA_S2S5S1S4vsS3S6S1S4_ECG = ft_math(cfg,GA_S2S5vsS1S4,GA_S3S6vsS1S4);


%from tutorial %works
cfg = [];
cfg.xlim = [0.472 0.484]; %0.4 0.5;0.4 0.468;0.4 0.452
%cfg.zlim = [-0.2 0];
cfg.layout = 'biosemi64.lay';
cfg.parameter = 'avg'; % the default 'avg' 'individual' is not present in the data
figure; ft_topoplotER(cfg,GA_S2S5vsS1S4); colorbar%GA_S2S5vsS1S4 %GA_S3S6vsS1S4

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
ft_singleplotER(cfg,GA_S1vsS4_ECG,GA_S2vsS5_ECG,GA_S3vsS6_ECG); %,GA_S4,GA_
