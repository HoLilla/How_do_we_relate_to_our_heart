%needs variables from study4_cpb_clean.m

close all
clear diffdat
clear mydat

mydat{1,1}      = GA_S1;
mydat{1,2}      = GA_S2;
mydat{1,3}      = GA_S3;
mydat{2,1}      = GA_S4;
mydat{2,2}      = GA_S5;
mydat{2,3}      = GA_S6;
mydat{1,1}.all  = iGA_S1.individual;
mydat{1,2}.all  = iGA_S2.individual; 
mydat{1,3}.all  = iGA_S3.individual;
mydat{2,1}.all  = iGA_S4.individual;
mydat{2,2}.all  = iGA_S5.individual;
mydat{2,3}.all  = iGA_S6.individual;
Rfrontal= {'AF4', 'AF8', 'F4', 'F6'}

% Parameters
TIME_WINDOW = [0.15 0.5]; % time in secs
CHANNELS    = Rfrontal;
BOOTSTRAPS  = 1000;
STRATEGY    = {'Attention','Perception','Regulation'};
CONG        = {'Congruent','Incongruent'};
CMAP        = cbrewer('qual','Set1', length(STRATEGY), 'cubic');

lay.label={'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';
'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';
'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';
'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';
'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';
'PO4';'O2';'COMNT';'SCALE'} 
% Generate indices to locate channels and time windows of interest
erp_time    = find(mydat{1}.time >= TIME_WINDOW(1) & mydat{1}.time <= TIME_WINDOW(2));
chans       = ismember(lay.label,CHANNELS);
figure('units','normalized','outerposition',[0 0 0.6 1.0]);
for a = 1:size(mydat,1) % COngruency
    
    for b = 1:size(mydat,2) % Strategy
        strat_dat = squeeze(mean(mydat{a,b}.all(:,chans,erp_time)));
        
        % Bootstrapping
        temp    = []; 
        error   = [];
        temp    = bootstrp(BOOTSTRAPS, @(x) mean(x), strat_dat);
        temp    = sort(temp);
        error(:,1)  = temp(975,:) - temp(500,:);
        error(:,2)  = temp(500,:) - temp(25,:); 
        
        % Subplot layout: Plotting Congruent data over Incongruent data
        % Change to subplot(size(mydat,1),2,2*a - 1) to include a
        % difference plot in column 2 (i.e cong minus incong with 95% CI)

        subplot(size(mydat,1),2,2*a-1)%2*a - 1         WORKS:subplot(size(mydat,1),1,a)
     

        % 95% CI Plots
        if b==1
        h{b} = boundedline(mydat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(3,:),'alpha','transparency',0.2);
        elseif b==2
        h{b} = boundedline(mydat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(2,:),'alpha','transparency',0.2);    
        elseif b==3
        h{b} = boundedline(mydat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(1,:),'alpha','transparency',0.2);
        end
        set(get(get(h{b}(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on
        
        % Plot of means
        if b==1
        p1=plot(mydat{1}.time(erp_time),mean(mydat{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(3,:));
        elseif b==2
        p2=plot(mydat{1}.time(erp_time),mean(mydat{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(2,:));
        elseif b==3
        p3=plot(mydat{1}.time(erp_time),mean(mydat{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(1,:));   
        end
        hold on
    end
   
    
    % Flash out the image with legends, titles and grids
    legend([p1 p2 p3],{'Attend','Feel','Regulate'}, 'location','northeast')
    %legend(axes_handle,'boxoff')
    p4=plot([mydat{1}.time(erp_time(1)) mydat{1}.time(erp_time(end))],[0 0],'LineStyle','--','Color','black');
   %middle line above
   set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
      p4b=plot([0.4 0.4],[-1.6 1.3],'LineStyle',':','Color','black');%hep from here line
    

   %middle line above
   set(get(get(p4b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
  
      hold on
   p5= plot([mydat{1}.time(erp_time(1)) mydat{1}.time(erp_time(end))],[-1.5 1],'LineStyle','none');
   set(get(get(p5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
  
%    %here i can specify the values 
    axis tight
    if a==1
     title({'AF4 AF8 F4 F6 - Congruent feedback'})
    elseif a==2
     title({'AF4 AF8 F4 F6 - Incongruent feedback'})
    
    end
     xlabel('Time (secs)')
     ylabel ('Amplitude (uV)')
     set(gca,'xtickLabel',{0:0.05:0.3})
end

hold on
clear diffdat

diffdat{1,1}     = GA_S2vsS1;
diffdat{2,1}     = GA_S3vsS1;
diffdat{3,1}     = GA_S3vsS2;
diffdat{4,1}     = GA_S5vsS4;
diffdat{5,1}     = GA_S6vsS4;
diffdat{6,1}     = GA_S6vsS5;
diffdat{1,1}.all = d_S2_S1.individual;
diffdat{2,1}.all = d_S3_S1.individual;
diffdat{3,1}.all = d_S3_S2.individual;
diffdat{4,1}.all = d_S5_S4.individual;
diffdat{5,1}.all = d_S6_S4.individual;
diffdat{6,1}.all = d_S6_S5.individual;

% Parameters
TIME_WINDOW = [0.15 0.5]; % time in secs
CHANNELS    = Rfrontal;
BOOTSTRAPS  = 1000;
STRATEGY    = {'Attention','Perception','Regulation'};
CONG        = {'Congruent','Incongruent'};
CMAP        = cbrewer('qual','Dark2', length(STRATEGY), 'cubic');
%CMAP        = cbrewer('seq','Greys', 6, 'cubic');

lay.label={'Fp1';'AF7';'AF3';'F1';'F3';'F5';'F7';'FT7';'FC5';'FC3';'FC1';
'C1';'C3';'C5';'T7';'TP7';'CP5';'CP3';'CP1';'P1';'P3';'P5';'P7';
'P9';'PO7';'PO3';'O1';'Iz';'Oz';'POz';'Pz';'CPz';'Fpz';'Fp2';'AF8';
'AF4';'AFz';'Fz';'F2';'F4';'F6';'F8';'FT8';'FC6';'FC4';'FC2';'FCz';
'Cz';'C2';'C4';'C6';'T8';'TP8';'CP6';'CP4';'CP2';'P2';'P4';'P6';'P8';'P10';'PO8';
'PO4';'O2';'COMNT';'SCALE'}
% Generate indices to locate channels and time windows of interest
erp_time    = find(diffdat{1}.time >= TIME_WINDOW(1) & diffdat{1}.time <= TIME_WINDOW(2));
chans       = ismember(lay.label,CHANNELS);

for b = 1:size(diffdat,1) % Strategy
       strat_dat = squeeze(mean(diffdat{b,1}.all(:,chans,erp_time)));
        
        % Bootstrapping
        temp    = []; 
        error   = [];
        temp    = bootstrp(BOOTSTRAPS, @(x) mean(x), strat_dat);
        temp    = sort(temp);
        error(:,1)  = temp(975,:) - temp(500,:);
        error(:,2)  = temp(500,:) - temp(25,:); 
        
        % Subplot layout: Plotting Congruent data over Incongruent data
        % Change to subplot(size(mydat,1),2,2*a - 1) to include a
        % difference plot in column 2 (i.e cong minus incong with 95% CI)
        subplot(6,size(mydat,1),2*b)%2*a - 1         WORKS:subplot(size(mydat,1),1,a)



        % 95% CI Plots
        h{b} = boundedline(diffdat{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(3,:),'alpha','transparency',0.2);
        set(get(get(h{b}(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on
        
        % Plot of means
        
        p1=plot(diffdat{1}.time(erp_time),mean(diffdat{b,1}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(3,:));
        hold on
    
   
    
    % Flash out the image with legends, titles and grids
    if b==1
    legend([p1],{'Feel-Attend (Congruent)'}, 'location','northwest')
    elseif b==2
    legend([p1],{'Regulate-Attend (Congruent)'}, 'location','northwest') 
    elseif b==3
    legend([p1],{'Regulate-Feel (Congruent)'}, 'location','northwest')
    elseif b==4
    legend([p1],{'Feel-Attend (Incongruent)'}, 'location','northwest')
    elseif b==5
    legend([p1],{'Regulate-Attend (Incongruent)'}, 'location','northwest')
    else b==6
    legend([p1],{'Regulate-Feel (Incongruent)'}, 'location','northwest')
    end
    p2=plot([diffdat{1}.time(erp_time(1)) diffdat{1}.time(erp_time(end))],[0 0],'LineStyle','--','Color','black');
   %middle line above
   set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
      if b==5
      p4b=plot([0.4 0.4],[-0.5 1.1],'LineStyle',':','Color','black');
      p4c=plot([0.468 0.468],[-0.5 1.1],'LineStyle',':','Color','black');

   set(get(get(p4b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
   set(get(get(p4c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      elseif b==4
      p4b=plot([0.4 0.4],[-0.5 1.1],'LineStyle',':','Color','black');
      p4c=plot([0.5 0.5],[-0.5 1.1],'LineStyle',':','Color','black');

   %middle line above
   set(get(get(p4b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
   set(get(get(p4c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      end
   hold on
   p3= plot([diffdat{1}.time(erp_time(1)) diffdat{1}.time(erp_time(end))],[-0.5 0.5],'LineStyle','none');
   set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
    
%    %here i can specify the values 
    axis tight
    
     xlabel('Time (secs)')
     ylabel ('Amplitude (uV)')
     set(gca,'xtickLabel',{0:0.05:0.3})
end


