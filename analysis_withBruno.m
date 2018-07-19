%% Analysis with Bruno
%%
%% path and toolboxes
analysisPath = '/home/raid2/nierula/Documents/OtherProjects/NoselfP3b/Analysis/EEG/';
% bbci toolbox
cd('/home/raid2/nierula/Documents/MATLAB/bbci_public-master/')
startup_bbci_toolbox()
% eeglab
addpath('/home/raid2/nierula/Documents/MATLAB/eeglab14_1_1b/')
% scripts
addpath(genpath([analysisPath 'Scripts']))

%% Preprocessing step: get epochs
%  ===================================== 
for experiment = 1:2  

    if experiment == 1
        myPath = [analysisPath 'ExpJune/'];
        nSubjects = [1:2 4:25];
    elseif experiment == 2
        myPath = [analysisPath 'ExpOctober/'];
        nSubjects = 26:35;
    end
    cd(myPath)

    for subject = nSubjects
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end

        %% load epo
        load cnt_raw cnt mrk

        %% make epochs for unsmoothed data
        ivEpo = [-200 800];
        epo = makeEpochs(cnt,mrk,ivEpo);

        save epo_raw epo
        clear epo

        %% smooth data
        [b,a] = butter(2,[0.1 10]/(cnt.fs/2));
        cnt = proc_filtfilt(cnt,b,a);        %% into 2 directions (back and forward)

        %% make epochs
        epo = makeEpochs(cnt,mrk,ivEpo);

        save epo_smoothed epo
    end
end

%% Preprocessing step: identification of artifacts
%  ================================================ 
clearvars -except analysisPath
for experiment = 1:2  

    if experiment == 1
        myPath = [analysisPath 'ExpJune/'];
        nSubjects  =  [1:2 4:25];
    elseif experiment == 2
        myPath = [analysisPath 'ExpOctober/'];
        nSubjects  =  26:35;
    end
    cd(myPath)
    
    for subject = nSubjects
        subject
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end
        
        %% load epo
        load epo_smoothed epo
                
        %% identify bad epochs by their amplitude being > 75 (on the absolute signal)
        for channel = 1:length(epo.clab)
            for trials = 1:size(epo.x,3)
                maxAmplitude(channel,trials) = max(abs(epo.x(:,channel,trials)));
            end
        end
        
        ampHigher75  =  maxAmplitude > 75;
        % if an epoch exceeds 75uV in one channel we interpolate this channel in this epoch only,
        % otherwise if this happens in more than one channel we remove the whole epoch from the analysis
        removedEpoch  =  [];
        interpolEpoch  =  [];
        interpolChannel  =  [];
        count1  =  0; count2  = 0;
        for epoch  =  1:size(ampHigher75, 2)
            if sum(ampHigher75(:, epoch))  ==  1 
                % interpolate epoch for this channel
                count1  =  count1 + 1;
                interpolEpoch(count1)  =  epoch;
                chan  =  find(ampHigher75(:, epoch)  ==  1);
                interpolChannel(count1)  =  chan; clear chan
            elseif sum(ampHigher75(:, epoch)) > 1
                % remove epoch from analysis
                count2  =  count2 + 1;
                removedEpoch(count2)  =  epoch;
            end
        end
        
        art.epoNot  =  removedEpoch;
        art.epoInterpol  =  interpolEpoch; 
        art.chanInterpol  =  interpolChannel;
    
        % save
        save artifacts art maxAmplitude
        
        clearvars -except myPath experiment subject analysisPath
        
    end
end

%% Preprocessing step: removal or interpolation of identified bad epochs
%  ====================================================================== 
clearvars -except analysisPath
for signal  =  1:2 % for smoothed (1) and raw (2) data
    
    for experiment = 1:2  

        if experiment == 1
            myPath = [analysisPath 'ExpJune/'];
            nSubjects  =  [1:2 4:25];
        elseif experiment == 2
            myPath = [analysisPath 'ExpOctober/'];
            nSubjects  =  26:35;
        end
        cd(myPath)

        for subject = nSubjects
            subject
            if subject<10
                cd([myPath 'S0' num2str(subject)])
            else
                cd([myPath 'S' num2str(subject)])
            end

            %% load data
            if signal  ==  1
                load epo_smoothed epo
            else 
                load epo_raw epo
            end
            load artifacts art

            %% interpolate epochs
            intEpo  =  art.epoInterpol;
            chanIdx  =  art.chanInterpol;

            for i  =  1:length(intEpo)
                epo1  =  proc_selectEpochs(epo, intEpo(i));
                epo1  =  interpolateChannels(epo1, chanIdx(i));
                epo.x(:, chanIdx(i), intEpo(i))  =  epo1.x(:,chanIdx(i));
                clear epo1
            end

            %% remove epochs
            remEpo  =  art.epoNot;
            epo  =  proc_selectEpochs(epo, 'not', remEpo);

            %% save 
            if signal  ==  1
                save epo_smoothed_clean75 epo
            else
                save epo_clean75 epo
            end
        end
    end
end


%% identify P3 range on grand average (on go-nogo difference curve)
%  ================================================================= 
clearvars -except analysisPath 
for experiment = 1:2  

    if experiment == 1
        myPath = [analysisPath 'ExpJune/'];
        nSubjects = [1:2 4:25];
    elseif experiment == 2
        myPath = [analysisPath 'ExpOctober/'];
        nSubjects = 26:35;
    end
    cd(myPath)

    for subject = nSubjects
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end

        %% load epo
        load epo_smoothed_clean75 epo

        %% take average and set to baseline
        ivBsline = [-200 0];
        epoav = proc_average(epo); clear epo
        epoav = proc_baseline(epoav,ivBsline);
    
        %% calculate time window around P3 peak for averaging 
        epoav = proc_selectChannels(epoav,'Pz');

        %1) average all three conditions
        epoav1 = proc_selectClasses(epoav,1:3);  %the 3 nogo conditions --> induce the P300
        epoav2 = proc_selectClasses(epoav,4:6); clear epoav  %the 3 go conditions
        epoav1.x = mean(epoav1.x,3);
        epoav2.x = mean(epoav2.x,3);
        % difference curve between nogo and go contitions
        epoav = epoav1;
        epoav.x = epoav1.x-epoav2.x; clear epoav1 epoav2
        epoav.y = 1; epoav.className = {'difference wave of mean nogo and mean go conditions'};
        
        %% merge to one variable for all subjects
        if subject == 1
            all = epoav;
            all.x = NaN(length(epoav.t),length(epoav.clab),35);
            all.y = 1;
            all.className = {'averaged epo of all subects'};
        end
        all.x(:,:,subject) = epoav.x;
        clearvars -except all subject experiment analysisPath myPath
    end
end

%% average participants and plot
cd(analysisPath)
load variables selectedSubj

allav = all; allav.x = [];
allav.x(:,1) = nanmean(all.x(:,1,selectedSubj),3);    %mean after removing nan values

% find p3 peaks and check in plot if they are correct
figure; plot(allav.t,allav.x(:,1,1)); hold on;

%2) get P3 peak latency (in difference curve)
x = allav.x;
[~,idx] = findpeaks(x(:,1));
[~,idx1] = max(allav.x(idx)); clear x
idx2 = idx(idx1);clear idx idx1
plot(allav.t(idx2),allav.x(idx2,1),'*')

%3) get N2 peak latency
x = allav.x*(-1);
[~,idx] = findpeaks(x(:,1));
idx1 = idx(find(allav.t(idx)<allav.t(idx2))); clear idx x
plot(allav.t(idx1(end)),allav.x(idx1(end),1),'*')

%4) calculate mid-latency (ml) between N2 and P3
ml = round((idx2-idx1(end))/2);

%5) define time window: P3-latency +/- ml
p3range = [round(allav.t(idx2-ml)) round(allav.t(idx2+ml))];

% save
save p3_range p3range


%% identify N1/P1 ranges in grand average curve over all subjects
%  =============================================================== 
clearvars -except analysisPath
for experiment = 1:2  

    if experiment == 1
        myPath = [analysisPath 'ExpJune/'];
        nSubjects = [1:2 4:25];
    elseif experiment == 2
        myPath = [analysisPath 'ExpOctober/'];
        nSubjects = 26:35;
    end
    cd(myPath)
    for subject = nSubjects
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end

        %% load epo
        load epo_clean75 epo

        %% take average and set to baseline
        ivBsline = [-200 0];
        epoav = proc_average(epo); clear epo
        epoav = proc_baseline(epoav,ivBsline);

        %% define time window around N1/P1 peaks for averaging 
        epoav = proc_selectChannels(epoav,{'PO3' 'PO4'});

        %% average over all conditions
        epoav.x = mean(epoav.x,3);

        %% merge to one variable for all subjects
        if subject == 1
            all = epoav;
            all.x = NaN(length(epoav.t),length(epoav.clab),35);
            all.y = 1;
            all.className = {'averaged epo of all subects'};
        end
        all.x(:,:,subject) = epoav.x;
        clearvars -except counter all subject myPath experiment analysisPath
    end
end

%% average participants and plot
cd(analysisPath)
load variables selectedSubj

allav = all; allav.x = [];
allav.x(:,1) = nanmean(all.x(:,1,selectedSubj),3);    %mean after removing nan values
allav.x(:,2) = nanmean(all.x(:,2,selectedSubj),3);    %mean after removing nan values

% find n1 and p1 peaks and check in plot if they are correct
figure, plot(allav.t,allav.x(:,1),'b'), hold on,
plot(allav.t,allav.x(:,2),'r'),legend(allav.clab),
xlabel('time [ms]'), ylabel('Amplitude [\muVolt]'),
title(['Experiment' num2str(experiment) ', N = ' num2str(size(all.x,3))]);

[~,maxIdx]  =  findpeaks(allav.x(:,1));
dataInv  =  1.01*max(allav.x(:,1)) - allav.x(:,1);
[~,minIdx]  =  findpeaks(dataInv);
idx1 = sort([maxIdx; minIdx]); clear maxIdx minIdx dataInv
peaks = find(allav.t(idx1)>100 & allav.t(idx1)<250);
if isempty(peaks) 
    % select peak with mouse
    % todo
    peaks  =  find(allav.t > 165);
    p1(1,1)  =  peaks(1);
    peaks  =  find(allav.t > 177);
    n1(1,1) = peaks(1); clear peaks
else
    p1(1,1) = idx1(peaks(1));
    n1(1,1) = idx1(peaks(2)); clear peaks
end

plot(allav.t([n1(1,1) p1(1,1)]),allav.x([n1(1,1) p1(1,1)],1),'b*')

[~,maxIdx]  =  findpeaks(allav.x(:,2));
dataInv  =  1.01*max(allav.x(:,2)) - allav.x(:,2);
[~,minIdx]  =  findpeaks(dataInv);
idx2 = sort([maxIdx; minIdx]); clear maxIdx minIdx dataInv
peaks = find(allav.t(idx2)>100 & allav.t(idx2)<250);
p1(1,2) = idx2(peaks(1));
n1(1,2) = idx2(peaks(2)); clear peaks
plot(allav.t([n1(1,2) p1(1,2)]),allav.x([n1(1,2) p1(1,2)],2),'r*')

%% mid-latency between n1 and p1
X  =  allav.t(n1) - allav.t(p1);
midLat  =  X/2; % time window: n1Max +/- midLat, p1Max +/- midLat,

%% P1 range
tpo3 = allav.t(p1(1,1));
rangePO3  =  [tpo3-midLat(1) tpo3+midLat(1)];
tpo4 = allav.t(p1(1,2));
rangePO4  =  [tpo4-midLat(2) tpo4+midLat(2)];
% save
save p1_range rangePO3 rangePO4
clear range* tpo*

%% N1 range
tpo3 = allav.t(n1(1,1));
rangePO3  =  [tpo3-midLat(1) tpo3+midLat(1)];
tpo4 = allav.t(n1(1,2));
rangePO4  =  [tpo4-midLat(2) tpo4+midLat(2)];
% save
save n1_range rangePO3 rangePO4 


%% get P3b values
% ================ 
clearvars -except analysisPath
load p3_range
for experiment = 1:2  

    if experiment == 1
        myPath = [analysisPath 'ExpJune/'];
        nSubjects = [1:2 4:25];
    elseif experiment == 2
        myPath = [analysisPath 'ExpOctober/'];
        nSubjects = 26:35;
    end
    cd(myPath)

    for subject = nSubjects
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end

        %% load epo
        load epo_smoothed_clean75 epo

        % take average and set to baseline
        ivBsline = [-200 0];
        epoav = proc_average(epo);
        epoav = proc_baseline(epoav,ivBsline);

        % save
        save epo_smoothed_clean75 epoav -append

        %% select channel
        epoav1 = proc_selectChannels(epoav,'Pz');

        %% select previously definde interval
        iv = p3range;%(subject,[2 3]);
        epoav1 = proc_selectIval(epoav1,iv);

        %% calculate mean amplitude in defined interval
        for condition = 1:6
            p3b(subject,condition) = mean(epoav1.x(:,1,condition));
        end
        clearvars -except p3b subject p3range myPath analysisPath
    end
end
cd(myPath)
cd ../
save p3 p3b

%% Plot grant avarage
clearvars -except analysisPath
cd(analysisPath)
load variables selectedSubj
load p3
counter = 0;
for experiment = 1:2  

    if experiment == 1
        myPath = [analysisPath 'ExpJune/'];
        nSubjects = selectedSubj(selectedSubj<26);
    elseif experiment == 2
        myPath = [analysisPath 'ExpOctober/'];
        nSubjects = selectedSubj(selectedSubj>25);
    end
    cd(myPath)
    for subject = nSubjects
        counter = counter+1;
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end
        load epo_smoothed_clean75 epoav
        epoav = proc_selectChannels(epoav,'Pz');
        if counter == 1
            epoav_all = epoav;
            epoav_all.x = [];
            epoav_all.title  =  'UnselfP3b all subjects';
        end
        epoav_all.x(:,:,:,counter)  =  epoav.x;
        epoav_all.N(counter,:)  =  epoav.N;
        clear epoav
    end
end
dat = epoav_all;
dat.x = [];
dat.x = mean(epoav_all.x,4);

figure;
plot(dat.t,dat.x(:,1,1),'k','LineWidth',2)
hold on;
plot(dat.t,dat.x(:,1,2),'r','LineWidth',2)
plot(dat.t,dat.x(:,1,3),'b','LineWidth',2)

plot(dat.t,dat.x(:,1,4),'--k','LineWidth',2)
plot(dat.t,dat.x(:,1,5),'--r','LineWidth',2)
plot(dat.t,dat.x(:,1,6),'--b','LineWidth',2)

title('Experiments 1 and 2'); axis([-200 800 -5 15]); legend(dat.className);
xlabel('time [ms]'); ylabel('P3b amplitude at Pz [\muV]');

%% get P1 and N1 values
% ====================== 
clearvars -except analysisPath
cd(analysisPath)
for erp = 1:2
    if erp == 1
        load p1_range
    else
        load n1_range
    end
    for experiment = 1:2  

        if experiment == 1
            myPath = [analysisPath 'ExpJune/'];
            nSubjects = [1:2 4:25];
        elseif experiment == 2
            myPath = [analysisPath 'ExpOctober/'];
            nSubjects = 26:35;
        end
        cd(myPath)

        for subject = nSubjects
            if subject<10
                cd([myPath 'S0' num2str(subject)])
            else
                cd([myPath 'S' num2str(subject)])
            end

            %% load epo
            load epo_clean75 epo

            % take average and set to baseline
            ivBsline = [-200 0];
            epoav = proc_average(epo);
            epoav = proc_baseline(epoav,ivBsline);
            
            % save
            save epo_clean75 epoav -append

            %% select channel
            epoav1 = proc_selectChannels(epoav,'PO3');
            epoav2 = proc_selectChannels(epoav,'PO4');

            %% select previously definde interval
            iv = rangePO3;%(subject,[2 3]);
            epoav1 = proc_selectIval(epoav1,iv); clear iv
            iv = rangePO4;%(subject,[2 3]);
            epoav2 = proc_selectIval(epoav2,iv); clear iv

            %% calculate mean amplitude in defined interval
            for condition = 1:6   % 1 = NoGo classic; 2 = NoGo syncMove; 3 = NoGo noMove; 
                                % 4 = Go classic; 5 = Go syncMove; 6 = Go noMove
                if ~ isempty(epoav1.x(:,:,condition))
                    po3(subject,condition) = mean(epoav1.x(:,:,condition));
                else
                    po3(subject,condition) = NaN;
                end
                if ~isempty(epoav2.x(:,:,condition))
                    po4(subject,condition) = mean(epoav2.x(:,:,condition));
                else
                    po4(subject,condition) = NaN;
                end
            end
            clearvars -except po3 po4 subject range myPath rangePO3 rangePO4 erp analysisPath
        end
    end
    cd(myPath)
    cd ../
    po3(3,:)  =  NaN; % to substitute 0 for missing subject #3
    po4(3,:)  =  NaN; % to substitute 0 for missing subject #3
    if erp == 1
        save p1 po3 po4
    else
        save n1 po3 po4
    end
    clear po3 po4
end