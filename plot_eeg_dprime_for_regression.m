%needs variables from study4_cpb_clean.m
%earlier filename: lh_sn_ERPplots_diffDiffAbsCombinded_regressionHEP.m 
clear diffdat

diffdat{1,1}      = GA_S1vsS4;
diffdat{1,2}      = GA_S2vsS5;
diffdat{2,1}      = GA_S1vsS4;
diffdat{2,2}      = GA_S3vsS6;
diffdat{1,1}.all  = d_S1_S4.individual;
diffdat{1,2}.all  = d_S2_S5.individual; 
diffdat{2,1}.all  = d_S1_S4.individual;
diffdat{2,2}.all  = d_S3_S6.individual;

% %6 rois sel et al
% Lfrontal=	{'AF7', 'AF3', 'F5', 'F3'}
% Lcentral= {'FC5', 'FC3', 'C5', 'C3'}
% Lparietal=	{'CP5', 'CP3', 'P5', 'P3'}
% 
% Rfrontal= {'AF4', 'AF8', 'F4', 'F6'}
% Rcentral= {'FC4', 'FC6', 'C4', 'C6'}
% Rparietal= {'CP4', 'CP6', 'P4', 'P6'}
% Parameters
TIME_WINDOW = [0.15 0.5]; % time in secs
CHANNELS    = Lcentral;
BOOTSTRAPS  = 1000;
STRATEGY    = {'Attention','Perception','Regulation'};
CMAP        = cbrewer('qual','Set1', length(STRATEGY), 'cubic');
%CMAP        = cbrewer('seq','Greys', 6, 'cubic');

%save diffdat diffdat
% Load necessary files

lay.label={'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';
'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';
'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';
'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';
'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';
'PO4';'O2';'COMNT';'SCALE'} 
% Generate indices to locate channels and time windows of interest
erp_time    = find(diffdat{1}.time >= TIME_WINDOW(1) & diffdat{1}.time <= TIME_WINDOW(2));
chans       = ismember(lay.label,CHANNELS);

figure('units','normalized','outerposition',[0 0 0.6 1.0]);
for a = 1:size(diffdat,1) % Strategy
    
    for b = 1:size(diffdat,2) % Congruence
        strat_dat = squeeze(mean(diffdat{a,b}.all(:,chans,erp_time)));
        
        % Bootstrapping
        temp    = []; 
        error   = [];
        temp    = bootstrp(BOOTSTRAPS, @(x) mean(x), strat_dat);
        temp    = sort(temp);
        error(:,1)  = temp(975,:) - temp(500,:);
        error(:,2)  = temp(500,:) - temp(25,:); 
        
        % Subplot layout: Plotting Congruent data over Incongruent data
        % Change to subplot(size(diffdat,1),2,2*a - 1) to include a
        % difference plot in column 2 (i.e cong minus incong with 95% CI)

        subplot(size(diffdat,1),2,2*a - 1)%2*a - 1         WORKS:subplot(size(diffdat,1),1,a)
     

        % 95% CI Plots
        if a==1 && b==1 
        h{b} = boundedline(diffdat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(3,:),'alpha','transparency',0.2);
        elseif a==1 && b==2 
        h{b} = boundedline(diffdat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(2,:),'alpha','transparency',0.2);
        elseif a==2 && b==1 
        h{b} = boundedline(diffdat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(3,:),'alpha','transparency',0.2);
        elseif a==2 && b==2 
        h{b} = boundedline(diffdat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(1,:),'alpha','transparency',0.2);
        end
        set(get(get(h{b}(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on
        
        % Plot of means
        if a==1 && b==1 
        p1=plot(diffdat{1}.time(erp_time),mean(diffdat{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(3,:));
        elseif a==1 && b==2
        p2=plot(diffdat{1}.time(erp_time),mean(diffdat{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(2,:));
        elseif a==2 && b==1
        p1=plot(diffdat{1}.time(erp_time),mean(diffdat{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(3,:));  
        elseif a==2 && b==2
        p2=plot(diffdat{1}.time(erp_time),mean(diffdat{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(1,:));            
        end
        hold on
    end
   
    
    % Flash out the image with legends, titles and grids
    p3=plot([diffdat{1}.time(erp_time(1)) diffdat{1}.time(erp_time(end))],[0 0],'LineStyle','--','Color','black');
   %middle line above
   set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
      p4b=plot([0.472 0.472],[-0.5 0.5],'LineStyle',':','Color','black');%hep from here line
      p4c=plot([0.484 0.484],[-0.5 0.5],'LineStyle',':','Color','black');%hep from here line

   %middle line above
   set(get(get(p4b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
   set(get(get(p4c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
      
      
   p4= plot([diffdat{1}.time(erp_time(1)) diffdat{1}.time(erp_time(end))],[-0.5 0.5],'LineStyle','none');
   set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
 
%    %here i can specify the values 
    axis tight
    if a==1
     title({'FC5 FC3 C5 C3 - Congruency difference'})
     legend([p1 p2],{'Attend(C-IC)','Feel(C-IC)'}, 'location','northeast')

    elseif a==2
     title({'FC5 FC3 C5 C3  - Congruency difference'})
     legend([p1 p2],{'Attend(C-IC)','Regulate(C-IC)'}, 'location','northeast')

    end
     xlabel('Time (secs)')
     ylabel ('Amplitude (uV)')
      set(gca,'xtickLabel',{0:0.05:0.3})
end

hold on
clear diffdiffdat

diffdiffdat{1,1}     = GA_S2S5vsS1S4;
diffdiffdat{2,1}     = GA_S3S6vsS1S4;
diffdiffdat{1,1}.all = d_S2S5_S1S4.individual;
diffdiffdat{2,1}.all = d_S3S6_S1S4.individual;

% Parameters
TIME_WINDOW = [0.15 0.5]; % time in secs
CHANNELS    = Rfrontal;
BOOTSTRAPS  = 1000;
STRATEGY    = {'Attention','Perception','Regulation'};
CMAP        = cbrewer('qual','Dark2', length(STRATEGY), 'cubic');
%CMAP        = cbrewer('seq','Greys', 6, 'cubic');

%save diffdat diffdat

lay.label={'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';
'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';
'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';
'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';
'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';
'PO4';'O2';'COMNT';'SCALE'} 
% Generate indices to locate channels and time windows of interest
erp_time    = find(diffdiffdat{1}.time >= TIME_WINDOW(1) & diffdiffdat{1}.time <= TIME_WINDOW(2));
chans       = ismember(lay.label,CHANNELS);

for a = 1:size(diffdiffdat,1) % Strategy
       strat_dat = squeeze(mean(diffdiffdat{a,1}.all(:,chans,erp_time)));
        
        % Bootstrapping
        temp    = []; 
        error   = [];
        temp    = bootstrp(BOOTSTRAPS, @(x) mean(x), strat_dat);
        temp    = sort(temp);
        error(:,1)  = temp(975,:) - temp(500,:);
        error(:,2)  = temp(500,:) - temp(25,:); 
        
        % Subplot layout: Plotting Congruent data over Incongruent data
        % Change to subplot(size(diffdat,1),2,2*a - 1) to include a
        % difference plot in column 2 (i.e cong minus incong with 95% CI)
        subplot(size(diffdat,1),2,2*a)%2*a - 1         WORKS:subplot(size(diffdat,1),1,a)



        % 95% CI Plots
        h{a} = boundedline(diffdiffdat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(3,:),'alpha','transparency',0.2);
        set(get(get(h{a}(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on
        
        % Plot of means
        
        p1=plot(diffdiffdat{1}.time(erp_time),mean(diffdiffdat{a,1}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(3,:));
        hold on
    
   
    
    % Flash out the image with legends, titles and grids
    
    p2=plot([diffdiffdat{1}.time(erp_time(1)) diffdiffdat{1}.time(erp_time(end))],[0 0],'LineStyle','--','Color','black');
   %middle line above
   set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
   p3= plot([diffdiffdat{1}.time(erp_time(1)) diffdiffdat{1}.time(erp_time(end))],[-0.5 0.5],'LineStyle','none');
   set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
    hold on
    p4b=plot([0.472 0.472],[-0.5 0.5],'LineStyle',':','Color','black');%hep from here line
      p4c=plot([0.484 0.484],[-0.5 0.5],'LineStyle',':','Color','black');%hep from here line
      

   %middle line above
   set(get(get(p4b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
   set(get(get(p4c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
      
      
%    %here i can specify the values 
    axis tight
    if a==1
     title({'Feel (C-IC) - Attend (C-IC)'})
     %legend([p1],{'Feel(C-IC)-Attend(C-IC)'}, 'location','northeast')
    elseif a==2
     title({'Regulate (C-IC) - Attend (C-IC)'})
     %legend([p1],{'Regulate(C-IC)-Attend(C-IC)'}, 'location','northeast')
    
    end
     xlabel('Time (secs)')
     ylabel ('Amplitude (uV)')
      set(gca,'xtickLabel',{0:0.05:0.3})
end
