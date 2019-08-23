close all
clear mydat_ECG
mydat_ECG{1,1}      = GA_S1_ECG;
mydat_ECG{1,2}      = GA_S2_ECG;
mydat_ECG{1,3}      = GA_S3_ECG;
mydat_ECG{2,1}      = GA_S4_ECG;
mydat_ECG{2,2}      = GA_S5_ECG;
mydat_ECG{2,3}      = GA_S6_ECG;
mydat_ECG{1,1}.all  = iGA_S1_ECG.individual;
mydat_ECG{1,2}.all  = iGA_S2_ECG.individual; 
mydat_ECG{1,3}.all  = iGA_S3_ECG.individual;
mydat_ECG{2,1}.all  = iGA_S4_ECG.individual;
mydat_ECG{2,2}.all  = iGA_S5_ECG.individual;
mydat_ECG{2,3}.all  = iGA_S6_ECG.individual;
ECG= {'EXG8'}

% Parameters
TIME_WINDOW = [0.15 0.5]; % time in secs
CHANNELS    = ECG;
BOOTSTRAPS  = 1000;
STRATEGY    = {'Attention','Perception','Regulation'};
CONG        = {'Congruent','Incongruent'};
CMAP        = cbrewer('qual','Set1', length(STRATEGY), 'cubic');

%save mydat mydat

lay.label={'EXG8'}
% Generate indices to locate channels and time windows of interest
erp_time    = find(mydat_ECG{1}.time >= TIME_WINDOW(1) & mydat_ECG{1}.time <= TIME_WINDOW(2));
chans       = ismember(lay.label,CHANNELS);
%figure('units','normalized','outerposition',[0 0 0.6 1.0]);
%figure('units','normalized','outerposition',[0 0 0.6 1.0]);
for a = 1:size(mydat_ECG,1) % COngruency
    
    for b = 1:size(mydat_ECG,2) % Strategy
        strat_dat = squeeze(mean(mydat_ECG{a,b}.all(:,chans,erp_time)));
        
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

        %subplot(size(mydat_ECG,1),2,2*a-1)%2*a - 1         WORKS:subplot(size(mydat,1),1,a)
     

        % 95% CI Plots
        if b==1
        h{b} = boundedline(mydat_ECG{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(3,:),'alpha','transparency',0.2);
        elseif b==2
        h{b} = boundedline(mydat_ECG{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(2,:),'alpha','transparency',0.2);    
        elseif b==3
        h{b} = boundedline(mydat_ECG{1}.time(erp_time),temp(500,:),error,'cmap',CMAP(1,:),'alpha','transparency',0.2);
        end
        set(get(get(h{b}(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        hold on
        
        % Plot of means
        if b==1
        if a==1
            p1=plot(mydat_ECG{1}.time(erp_time),mean(mydat_ECG{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(3,:));
        else
            p1b=plot(mydat_ECG{1}.time(erp_time),mean(mydat_ECG{a,b}.avg(chans,erp_time),1),'LineWidth',2,'LineStyle','--','Color',CMAP(3,:));
        end
        elseif b==2
            if a==1
        p2=plot(mydat_ECG{1}.time(erp_time),mean(mydat_ECG{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(2,:));
            else
         p2b=plot(mydat_ECG{1}.time(erp_time),mean(mydat_ECG{a,b}.avg(chans,erp_time),1),'LineWidth',2,'LineStyle','--','Color',CMAP(2,:)); 
            end
        elseif b==3
            if a==1
        p3=plot(mydat_ECG{1}.time(erp_time),mean(mydat_ECG{a,b}.avg(chans,erp_time),1),'LineWidth',1.5,'Color',CMAP(1,:));
            else
        p3b=plot(mydat_ECG{1}.time(erp_time),mean(mydat_ECG{a,b}.avg(chans,erp_time),1),'LineWidth',2,'LineStyle','--','Color',CMAP(1,:));
            end
        end
        hold on
   
    end
     
   
    
    
    p4=plot([mydat_ECG{1}.time(erp_time(1)) mydat_ECG{1}.time(erp_time(end))],[0 0],'LineStyle','--','Color','black');
   %middle line above
   set(get(get(p4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
      p4b=plot([0.4 0.4],[-2.8 12.1],'LineStyle',':','Color','black');%hep from here line
    % p4c=plot([0.5 0.5],[-1.6 1.3],'LineStyle',':','Color','black');%hep from here line

   %middle line above
   set(get(get(p4b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
   %set(get(get(p4c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
      hold on
   p5= plot([mydat_ECG{1}.time(erp_time(1)) mydat_ECG{1}.time(erp_time(end))],[-1.5 1],'LineStyle','none');
   set(get(get(p5,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');%turns it off
  
%      yerp_max = max([mydat{1}.time(erp_time(1)),max(mydat{1}.time(erp_time(1)))]);
%      yerp_min = min([mydat{1}.time(erp_time(1)),max(mydat{1}.time(erp_time(1)))]);
%    plot([0.1 0.5],[0.1 0.5],'LineStyle','--','Color','black');
%    %here i can specify the values 
    axis tight
    
     title({'ECG signal'})
    xlabel('Time (secs)')
     ylabel ('Amplitude (uV)')
     set(gca,'xtickLabel',{0:0.05:0.3})
end
hold on
%l1=legend([p1 p2 p3 p1b p2b p3b],{'Attend (Congruent)','Feel (Congruent)','Regulate (Congruent)','Attend (Incongruent)','Feel (Incongruent)','Regulate (Incongruent)'},'Location','south','NumColumns',2)
