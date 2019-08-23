clear all
close all

subjid=input('subject id:','s');

datapath=[pwd '\DATA\'];%generates new folder in current dir
age=input('ptp age: ' ,'s');
gender=input('ppt gender: ' ,'s');
device=input('device id:','s');

filename=[subjid 'Study4_LH'];
if exist ([datapath filename '.mat'],'file')
    filename =[filename (num2str(round(sum(fix(clock)))))]
end

dem(1).filename=filename;
dem(1).age=age;
dem(1).gender=gender;
%% Prepare asynch ryhtmns database
samples=[959,953,913,872,858,850,835,818,812,799,791,786,772,771,763,761,...
    734,725,697,695,691,684,683,676,662,640,639,634,594,1220,1129,974,924,880,837,830,...
    826,813,808,783,766,762,758,753,749,748,737,735,730,727,721,709,705,700,679,661,...
    638,628,624,603,598,555,568,881,807,831,508,936,894,1012,1176,1162];
for i=1:length(samples) %make it to same size - beat matrix
    clear ibi
    name=[num2str(samples(i))]; %loading matlab file with name
    load (name);
    IBI(i,:)=ibi(1:800); %makes i same length as takes the first 100 % I need longer than 100
end

%% define the display settings and initialize cogent(toolbox to present stimuli)
global cogent;
global ratingCounting;
global ratingDetect;

config_keyboard
config_display(1,1,[0 0 0], [1 1 1], 'Arial',20,12,0)
start_cogent

%% define variables and prepare stuff in the background
HB_recognition=[];
trialN=52;%%40 for final behavioural pilot+6 trials/subcondition to have >=250 hbs/condition for people around 60bpm
totalN=3*trialN; %for three conditions (count, feel, regulate)
nhb= 3; % HR estimates will consist the average IBI of the preceeding nhb heartbeats.
%The wider the smoother are the changes
nhbbase=10;%for baseline

rng('Shuffle') %randomises the seed, I only need it once in the script
%randomising the asynch and synch trials 1=synch, 2=asynch

Contingency_mtx=repmat([ones(1,trialN/2) ones(1,trialN/2)*2],3);
Contingency_mtx=Contingency_mtx(1,:);
slow_fast_mtx=repmat([ones(1,trialN/2)*0 ones(1,trialN/4) ones(1,trialN/4)*2],3);%1=fast 2=slow 0=nochange
slow_fast_mtx=slow_fast_mtx(1,:);

%randomising the startegy trials 1=count, 2=feel, 3=regulate
Strategy_mtx=[ones(1,trialN) ones(1,trialN)*2, ones(1,trialN)*3 ];
Condition_mtx=[slow_fast_mtx;Contingency_mtx;Strategy_mtx;];
index=randperm(length(Condition_mtx));

gf_r=[1 1 1 1 1 0 0 0 0 0 ];
gf_r=repmat(gf_r,1,10);
for t=1:length(Strategy_mtx)
    phdata(t).Condition=Condition_mtx(:,index(t));
    phdata(t).slowfast_order=phdata(t).Condition(1,:);
    phdata(t).Contingency_cond=phdata(t).Condition(2,:);
    phdata(t).Strategy_cond=phdata(t).Condition(3,:);
    phdata(t).gf=gf_r(randperm(length(gf_r)));    
end

%% Create background Thermometer to be called each time (Sprite 1)

cgmakesprite(1,640,480,0,0,0) %create a diff sheet index numbe 1 here, size,
cgsetsprite(1) %i will use this-write on it
cgpencol(1,1,1)
cgpencol(1,1,1)
cgpenwid(2)
cgdraw(-12,-102,-12,102)%this draws lines in positions, this give the structure to the TM
cgdraw(12,-102,12,102)
cgdraw(-12,-102,12,-102)
cgdraw(-12,102,12,102)

cgdraw(-20,100,-15,100) %ticks on the termometer
cgdraw(-20,75,-15,75)
cgdraw(-20,50,-15,50)
cgdraw(-20,25,-15,25)
cgdraw(-40,0,-15,0)
cgdraw(-20,-100,-15,-100)
cgdraw(-20,-75,-15,-75)
cgdraw(-20,-50,-15,-50)
cgdraw(-20,-25,-15,-25)

cgsetsprite(0)%index zero>>blank

%% This is to inatialize the communication protocol with the parallel port
% ioObj = io64;
% status = io64(ioObj);
% io64(ioObj,57336,0); %Change to to configure ports, each PC is different

if device=='g'    
    s=daq.createSession('ni');
    s.addDigitalChannel('Dev1','Port1/Line0','InputOnly');
    s.addDigitalChannel('Dev1','Port0/Line6','OutputOnly');%Output channel 1
    s.addDigitalChannel('Dev1','Port0/Line4','OutputOnly'); %Output channel2
    s.outputSingleScan([0 0]); %needs equal numbers of coloumns as outputchannels    
else
    s=daq.createSession('ni');
    s.addDigitalChannel('Dev2','Port1/Line0','InputOnly');
    s.addDigitalChannel('Dev2','Port0/Line6','OutputOnly');%Output channel 1
    s.addDigitalChannel('Dev2','Port0/Line4','OutputOnly'); %Output channel2
    s.outputSingleScan([0 0]); %needs equal numbers of coloumns as outputchannels
end

ioObj = io64;
status = io64(ioObj);
io64(ioObj,53504,0);
io64(ioObj,53504,128); wait(10); io64(ioObj,53504,0);

%% Instructions
cgpencol(1,1,1)
cgfont('Helvetica',20)
cgtext('Here is a short recap of the instructions before we start the real experiment',0,100)
cgpencol(0,0.8,0)
cgtext('COUNT:',0,80)
cgpencol(1,1,1)
cgtext('count the GREEN pulses',0,60)
cgpencol(0,0.5,1)
cgtext('FEEL:',0,30)
cgpencol(1,1,1)
cgtext('do you FEEL a heartbeat when a BLUE pulse is present?',0,10)
cgpencol(1,0,0)
cgtext('REGULATE:',0,-20)
cgpencol(1,1,1)
cgtext('bring DOWN the bar by relaxing',0,-40)

cgflip(0,0,0)
clearkeys;
readkeys;
[key ktime pressed]=waitkeydown(inf, 71); % wait infinite time for a spacebar press

%measure baseline
cgpencol(1,1,1)
cgfont('Helvetica',20)
cgtext('Please wait a few seconds.',0,0)
cgtext('Please remember to be very still',0,-50)
cgflip(0,0,0)

%% this will estimate the first HB
% it works like this: heartbeats are communicated from the Powerlab to Matlab
%through a parallel port
% this script reads the state of the input of the parallel port to know
% if it is currently "on" (heartbeat (R-wave) detected) or "off" (between
% hearbeats (R-waves)).
% if it is "on" the port will be active (binary 1)and have a value of 144
%(this changes from PC to PC and parellel card to parallel card)
% if it is "off" the port will be inactive (binary 0)a value of 128
%(this changes from PC to PC and parellel card to parallel card)
% remeber that the port will be active for the entire legnth of the pulse
% defined in powerlab (e.g. 20ms)

heart=0; % the input parallel port state %128 no HB
while heart==0  %% it will continuosly check the port state until it is not longer 128 (= binary 0). when this happens an hearbeat was detected
    heart=inputSingleScan(s); % check input port state and assign it to "heart"
end
t0=time; % an hearbeat was detected. get timestamp for this heartbeat
s.outputSingleScan([1 0]); wait(10); s.outputSingleScan([0 0]); %sends trigger - visualises trigger % send a pulse to powerlab just to visualize it (this is not necessary)

while heart~=0 %% wait for the pulse to finnish 144 is hb %%it will continuosly check the port state until it is not longer 144 (binary 1). when this happens the hearbeat pulse has finished
    heart=inputSingleScan(s);
end

%% measure thermomether for the first time without any signal to generate ibi
%series for non-contingent trials and first baseline
for b=1:nhbbase %% the variable "b" will track the number of heartbeats
    while heart==0
        heart=inputSingleScan(s); %% wait for an heartbeat
    end
    tHB(b)=time-t0; % calculate and save in the tHB variable the interbeat interval of heartbeat(b)
    t0=time; %save time of last heartbeat
    s.outputSingleScan([1 0]); wait(40); s.outputSingleScan([0 0]); %sends trigger - visualises trigger % send a pulse to powerlab just to visualize it (this is not necessary)
    
    while heart~=0
        heart=inputSingleScan(s);
    end
    
    meanIBI(b)=mean(tHB(1:b));  %estimate average interbeat interval for the previous n heartbeats
    meanHR(b)=60000/meanIBI(b); %estimate average heartrate for the previous n heartbeats
    previous_trial_IBI(b)=meanIBI(b); % for first time
    
end

    phdata(1).meanIBI(1)=meanIBI(b);
    phdata(1).level(1)=NaN;

%% main part
%%
self_reference_IBI= phdata(1).meanIBI(1);

fmatch=self_reference_IBI-samples; %samples are filenames and average of that ibi series
dif=sort(abs(fmatch));
MIN=dif(1);
x=find(abs(fmatch)==MIN);
ref=samples(x(1));

dem.OS_sample=ref;
OS_IBI=IBI(x,:);

pdif=((self_reference_IBI/ref)*100)-100; %percentage difference bw ref+other HR and my target
OSIbi=OS_IBI*(((pdif)/100)+1); %adjustment for being a diffferent heart
OSIbi=repmat(OSIbi,1,3); %looped 3times



start_NC_heart=20;
phdata(1).start_NC_heart=start_NC_heart;
%%
for trial_number = 1:totalN
    Strategy=phdata(trial_number).Strategy_cond
    Contingency=phdata(trial_number).Contingency_cond
    
    if trial_number~=1
        phdata(trial_number).start_NC_heart=phdata(trial_number-1).start_NC_heart+nHBnc;       
    end
    
    if Strategy==1
        cgpencol(0,0.8,0) %text is presented green
        cgfont('Helvetica',60)
        cgtext('COUNT!',0,0)
        cgflip(0,0,0)
        wait(1000)
        clearkeys;
        
    elseif Strategy==2
        cgpencol(0,0.5,1) %text is presented blue
        cgfont('Helvetica',60)
        cgtext('FEEL!',0,0)
        cgflip(0,0,0)
        wait(1000)
        clearkeys;
        
    elseif Strategy==3
        cgpencol(1,0,0) %text is presented red
        cgfont('Helvetica',60)
        cgtext('REGULATE!',0,0)
        cgflip(0,0,0)
        wait(1000)
        clearkeys;
    end
      
    
    %% sets the bar in the back ground and the scales for visualisin BF
    
    if trial_number==1
        previous_NC_IBI=mean(OSIbi(1:start_NC_heart));
    else
        previous_NC_IBI=mean(OSIbi(start_NC_heart-15:start_NC_heart)); %%baseline for non-contigent heart is always based on un-changed NC
    end
    NC_adjustment=phdata(trial_number).slowfast_order;
    if NC_adjustment==1
       phdata(trial_number).NC_heart=OSIbi(start_NC_heart:length(OSIbi))*0.85; %target average NC ibi faster - it is already looped so should be long enough
    phdata(trial_number).NC_heart_smooth=OSIbi(start_NC_heart-nhb+1:length(OSIbi))*0.85;
    elseif NC_adjustment==2
        phdata(trial_number).NC_heart=OSIbi(start_NC_heart:length(OSIbi))*1.15; %target average NC ibi slower
        phdata(trial_number).NC_heart_smooth=OSIbi(start_NC_heart-nhb+1:length(OSIbi))*1.15; %target average NC ibi slower
    elseif NC_adjustment==0
        phdata(trial_number).NC_heart=OSIbi(start_NC_heart:length(OSIbi));%no change
        phdata(trial_number).NC_heart_smooth=OSIbi(start_NC_heart-nhb+1:length(OSIbi));
    end
    
    for t=1:length(phdata(trial_number).NC_heart_smooth)-nhb
        phdata(trial_number).NC_heart_smooth(1,t)=mean(phdata(trial_number).NC_heart_smooth(t:t+nhb-1)); 
    end
    
  
    
    
    if Contingency==1
       phdata(trial_number).baseline=round(mean(previous_trial_IBI));  %estimate baseline, i.e. participants heart-rate in the beggining of the trial
    elseif Contingency==2
      %  phdata(trial_number).baseline=round(mean(previous_NC_IBI));  %estimate baseline, i.e. participants heart-rate in the beggining of the trial
         phdata(trial_number).baseline=phdata(trial_number).NC_heart(1);
    end
    scalemax=phdata(trial_number).baseline+(phdata(trial_number).baseline/2); %steps should be scaled to IBI
    scalemin=phdata(trial_number).baseline-(phdata(trial_number).baseline/4);
    steps=200/(phdata(trial_number).baseline);
    
    startscale=0; %% this is the initial value for the scale (min=0; max=200)
    
    cgdrawsprite(1,0,0)
    cgpencol(1,0,0) %RED COLOUR
    cgpenwid(20)
    cgdraw(0,-90,0,-10);
    cgflip(0,0,0)
    
    startTime=time;
    clearkeys;
    nHBnc=1;
    nHBc=1;
    k=0; %this variable will code for the space barpress (see below)
    
    heart=0; 
    while heart==0  %% it will continuosly check the port state until it is not longer 128 (= binary 0). when this happens an hearbeat was detected
        heart=inputSingleScan(s); % check input port state and assign it to "heart"
    end
    t0=time; % an hearbeat was detected. get timestamp for this heartbeat
    s.outputSingleScan([1 0]); wait(2); s.outputSingleScan([0 0]);
    
    while heart~=0
        heart=inputSingleScan(s);
    end
    
    tcont2=time;
    while k<10        
        contin=1; %continue
        while contin==1
            %% it will continuosly check the port state until it is not longer 0.
            if Contingency==2
                if time-tcont2>phdata(trial_number).NC_heart(nHBnc)
                    phdata(trial_number).HBnctime(nHBnc)=time;
                    nHBnc=nHBnc+1;
                    tcont2=time;
                    contin=0;
                end
            end
            heart=inputSingleScan(s); % check input port state and assign it to "heart"
            if heart~=0 && time-t0>100
                b=b+1; % increase number of detected heartbeats by one                               
                tHB(b)=time-t0;% calculate and save in the tHB variable the interbeat interval of heartbeat(b)
                phdata(trial_number).HBctime(nHBc)=time; 
                phdata(trial_number).HBcIBI(nHBc)=time-t0;                 
                nHBc=nHBc+1;
                t0=time;
                s.outputSingleScan([1 0]); wait(2); s.outputSingleScan([0 0]); ;% send a pulse to powerlab just to visualize it (this is not necessary)
                
                if Contingency==1 && Strategy==1
                    io64(ioObj,53504,1); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Contingent+count
                elseif Contingency==1 && Strategy==2
                    io64(ioObj,53504,2); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Contingent+feel
                elseif Contingency==1 && Strategy==3
                    io64(ioObj,53504,4); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Contingent+regulate
                elseif Contingency==2 && Strategy==1
                    io64(ioObj,53504,8); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Non-contingent+count
                elseif Contingency==2 && Strategy==2
                    io64(ioObj,53504,16); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Non-contingent+feel
                elseif Contingency==2 && Strategy==3
                    io64(ioObj,53504,32); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Non-contingent+regulate
                end
                if Contingency==1
                    contin=0;
                end
                
            end
        end
      
        meanIBI(b)=mean(tHB(b-nhb+1:b)); % average interbeat interval of the previous nhb heartbeats(nhb is defined at the top o fthe script)
        meanHR(b)=60000/meanIBI(b); % average heart rate of the previous nhb heartbeats(nhb is defined at the top o fthe script)
        phdata(trial_number).meanIBI(nHBc)=meanIBI(b);
        %in labchart, but this way we do not have to smooth it before feeding it in the thermometer level changes
        %so we don't need to transform tHB to meanIBI based on n=nhb consecutive heartbeats
        %current_trial_beat=current_trial_beat+1;
              
        if Contingency==1 %Contingent BF
            level=(round((meanIBI(b)-phdata(trial_number).baseline)*steps)*(-1)); %the scale is inverted (i.e, *(-1)) such that longer ibi= lower levels            
             phdata(trial_number).level(nHBc)=level;
        elseif Contingency==2 %Non-contingent          
            level=(round(phdata(trial_number).NC_heart_smooth(nHBnc)-phdata(trial_number).baseline)*steps)*(-1); %the scale is inverted (i.e, *(-1)) such that longer ibi= lower levels
            phdata(trial_number).level(nHBnc)=level;
        end
        
        %% visualises the pulses of different colour
        %I think this needs to be conditional/linked to the condition
        if Contingency==1
            pulse= phdata(trial_number).gf(nHBc);
        elseif Contingency==2
            pulse= phdata(trial_number).gf(nHBnc);
        end
        cgdrawsprite(1,0,0)
        if  pulse==0
            cgpencol(1,0.5,0)%orange
        elseif pulse==1
            if Strategy==1
                cgpencol(0,1,0)%green
            elseif Strategy==2
                cgpencol(0,0.5,1)%blue
            elseif Strategy==3
                cgpencol(0.9,0.9,0.9)%white
            end
        end
        cgpenwid(20)
        cgdraw(0,-90,0,startscale+level-10);
        cgflip(0,0,0)
        twait=time;
        while time-twait<150
            heart=inputSingleScan(s);
            if heart~=0 && time-t0>100
                b=b+1; % increase number of detected heartbeats by one
                tHB(b)=time-t0;% calculate and save in the tHB variable the interbeat interval of heartbeat(b)
                phdata(trial_number).HBctime(nHBc)=time; 
                phdata(trial_number).HBcIBI(nHBc)=time-t0;                 
                nHBc=nHBc+1;
                t0=time;
                s.outputSingleScan([1 0]); wait(2); s.outputSingleScan([0 0]); ;% send a pulse to powerlab just to visualize it (this is not necessary)
                
                if Contingency==1 && Strategy==1
                    io64(ioObj,53504,1); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Contingent+count
                elseif Contingency==1 && Strategy==2
                    io64(ioObj,53504,2); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Contingent+feel
                elseif Contingency==1 && Strategy==3
                    io64(ioObj,53504,4); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Contingent+regulate
                elseif Contingency==2 && Strategy==1
                    io64(ioObj,53504,8); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Non-contingent+count
                elseif Contingency==2 && Strategy==2
                    io64(ioObj,53504,16); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Non-contingent+feel
                elseif Contingency==2 && Strategy==3
                    io64(ioObj,53504,32); wait(2); io64(ioObj,53504,0);% send a pulse to actiview when Non-contingent+regulate
                end
            end
        end
        
        cgdrawsprite(1,0,0)
        cgpencol(1,0,0) %red
        cgpenwid(20)
        cgdraw(0,-90,0,startscale+level-10);
        cgflip(0,0,0)
          
        meanIBI(b)=mean(tHB(b-nhb+1:b)); % average interbeat interval of the previous nhb heartbeats(nhb is defined at the top o fthe script)
        meanHR(b)=60000/meanIBI(b); % average heart rate of the previous nhb heartbeats(nhb is defined at the top o fthe script)
        phdata(trial_number).meanIBI(nHBc)=meanIBI(b);
        
        readkeys;
        [key ktime pressed] = getkeydown; %check if any key was pressed
        clearkeys;
        
        if time-startTime>=10*1000 %10 seconds
            k=100 ;
        end
        if key==71  % if the spacebar was pressed end the loop
            k=100 ;
        end
    end
    %%%%%%%%%%%%%% end of trial %%%%%%%%%%%%%%
    
    s.outputSingleScan([1 0]); wait(100); s.outputSingleScan([0 0]);% thicker trigger to labchart as it ends
    %% behavioural response part   
    
    if Contingency==1
        phdata(trial_number).gf=phdata(trial_number).gf(1:nHBc-1);
        
    else
        phdata(trial_number).gf=phdata(trial_number).gf(1:nHBnc-1);
        start_NC_heart=start_NC_heart+nHBnc;  
        phdata(trial_number).HBncIBI=diff(phdata(trial_number).HBnctime);
    end   
    phdata(trial_number).performance_baseline=round(mean(previous_trial_IBI));
    previous_trial_IBI=mean(phdata(trial_number).HBcIBI);
           
    if Strategy==1
        take_ratings_count
        count_rating(1,1) = ratingCounting;
        phdata(trial_number).selfresponse=count_rating(1,1);
        ratingCounting=0;
        clearpict;
        %task_performance= greens
        greens=sum(phdata(trial_number).gf);
        if phdata(trial_number).selfresponse>=greens
            phdata(trial_number).performance=1-(phdata(trial_number).selfresponse/greens-1);
        else
            phdata(trial_number).performance=(phdata(trial_number).selfresponse/greens);
        end
        phdata(trial_number).pulse_number=greens;
        
    elseif Strategy==2
        take_ratings_feel2
        rating(1,1) = ratingDetect;
       phdata(trial_number).selfresponse=rating(1,1);
        ratingCounting=0;
        clearpict;
        blues=sum(phdata(trial_number).gf);       
        %task_performance
        if Contingency==1
            if(phdata(trial_number).selfresponse>5 && phdata(trial_number).selfresponse<45)
                phdata(trial_number).performance=(blues/2)/blues;
            elseif phdata(trial_number).selfresponse<5
                phdata(trial_number).performance=0;
            elseif phdata(trial_number).selfresponse>45
                phdata(trial_number).performance=1;
            end
        elseif Contingency==2
            if (phdata(trial_number).selfresponse>5 && phdata(trial_number).selfresponse<45)
                phdata(trial_number).performance=(blues/2)/blues;
            elseif phdata(trial_number).selfresponse<5
                phdata(trial_number).performance=1;
            elseif phdata(trial_number).selfresponse>45
                phdata(trial_number).performance=0;
            end
        end
        phdata(trial_number).pulse_number=blues;
        
    elseif Strategy==3
        take_ratings_regulate
        rating(1,1) = ratingDetect;
        phdata(trial_number).selfresponse=rating(1,1);
        ratingCounting=0;
        clearpict;
        whites=sum(phdata(trial_number).gf);
        %task performance decrease in HR from baseline
        phdata(trial_number).performance= 1-(mean(phdata(trial_number).HBcIBI)/phdata(trial_number).performance_baseline);
        phdata(trial_number).pulse_number=whites;
    end
    
    cgpencol(1,1,1) %text is presented in white
    cgfont('Helvetica',20)
    cgtext('Whose heart was the feedback representing?',0,0)
    cgtext('Not my heart                                         My heart',0,-60)
    cgflip(0,0,0)
    clearkeys;
    readkeys;
    [key ktime pressed]=waitkeydown(inf); % wait infinite time for any key to be pressed
    
    if key==4  % if any buton on the right was pressed on responsebox
        phdata(trial_number).HB_recognition=1; %1 is self
    elseif key==6
        phdata(trial_number).HB_recognition=1;
    elseif key==19 % if any buton on the left was pressed
        phdata(trial_number).HB_recognition=2; %2 is other
    elseif key==1 % if any buton on the left was pressed
        phdata(trial_number).HB_recognition=2;
    end
    
    take_ratings_Detection
    rating(1,1) = ratingDetect;
    phdata(trial_number).HB_confidence=rating(1,1);
    ratingDetect=0;
    clearpict;
    
    save([datapath filename],'phdata','dem');
    
    if trial_number == round(totalN/2)% to have break at half way through
        cgpencol(1,1,1)
        cgfont('Helvetica',20)
        cgtext('We are half way through the study.',0,0)
        cgtext('Please let the experimenter know you are ready to continue.',0,-20)
        cgflip(0,0,0)
        
        clearkeys;
        readkeys;
        [key ktime pressed]=waitkeydown(inf, 71); % wait infinite time for a spacebar press
        clearkeys;
    end
    
    cgfont('Helvetica',60)
    cgtext('+',0,0)
    cgflip(0,0,0)
    random_jitter=[750, 1250, 1750]
    index=randperm(length(random_jitter));
    wait(random_jitter(index))%jitter between random trials
end

for i=1:length(phdata) %1:correct 0:incorrect
    if phdata(i).Contingency_cond==phdata(i).HB_recognition
        phdata(i).correct_response(i)=1;
    else
        phdata(i).correct_response(i)=0;
    end
end

cgpencol(1,1,1)
cgfont('Helvetica',20)
cgtext('Thank you for taking part in the study!',0,0)
cgflip(0,0,0)
wait(3000)
%%
stop_cogent
close all

save([datapath filename],'phdata','dem');