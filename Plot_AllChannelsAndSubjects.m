path  =  '/home/raid2/nierula/Documents/OtherProjects/NoselfP3b/Analysis/EEG';
cd('/home/raid2/nierula/Documents/MATLAB/bbci_public-master/')
startup_bbci_toolbox()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ============  1) get average signal for included subjects  =============== 
counter = 0;
for experiment = 1:2  

    if experiment == 1
        myPath = [path 'ExpJune/'];
        nSubjects = [1 4 5 7 9:21 23 25]; % subjects 2, 3, 6, 8, 22, and 24 are already excluded
    elseif experiment == 2
        myPath = [path 'ExpOctober/'];
        nSubjects = [27:29 31:35]; % subjects 26, and 30 are already excluded
    end
    cd(myPath)

    for subject = nSubjects
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end
        counter = counter+1;
        
        % load epo
        load epo_smoothed_clean75 epoav
        % check if channel number is correct
        if size(epoav.x,2)<13 
            newEpoav = epoav;
            newEpoav.clab = []; newEpoav.x = [];
            lp = 1;
            for i = 1:13
                if size(epoav.clab,2)>=  i
                    if ~strcmp(epoav_all.clab{i},epoav.clab{lp})
                        newEpoav.clab{i} = epoav_all.clab{i};
                        newEpoav.x(:,i,:) = NaN;
                    else
                        newEpoav.clab{i} = epoav.clab{lp};
                        newEpoav.x(:,i,:) = epoav.x(:,lp,:);
                        lp = lp+1;
                    end
                else
                    newEpoav.clab{i} = epoav_all.clab{i};
                    newEpoav.x(:,i,:) = NaN;
                end
            end
            epoav = newEpoav;  
        end
    
        if counter == 1
            epoav_all = epoav;
            epoav_all.x = [];
        end
        epoav_all.x(:,:,:,counter) = epoav.x;
        clear epoav
    end
end

dat = epoav_all;
epoav.x = [];
epoav.x = mean(epoav_all.x,4);

cd(path)
save epoav_includedSubjects epoav

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ============  2) Plot ERP  =================================== 
clearvars -except path
cd(path)
load epoav_includedSubjects epoav

% all channels with subplots
gridPos =  [4 9 11 13 15 17 19 21 23 25 27 31 33];
chan =  epoav.clab;

% plot all channels
figure;
% xRange = input('Insert the x-Range in []');
% ax1 = [-200 800 xRange(1) xRange(2)];  % needs to be changed: ax1 = [xMin xMax yMin yMax]
ax1 = [-200 600 -5 15];  % needs to be changed: ax1 = [xMin xMax yMin yMax]

for i = 1:length(gridPos)
    subplot(5, 7, gridPos(i));hold on;
    plot(epoav.t,epoav.x(:, i, 1),'k') %% Critical: non-target
    hold on;
    %grid;
    plot(epoav.t,epoav.x(:, i, 2),'r') %% Critical: target
    plot(epoav.t,epoav.x(:, i, 3),'b') %% Control: non-target
    plot(epoav.t,epoav.x(:, i, 4),':k') %% Control: target
    plot(epoav.t,epoav.x(:, i, 5),':r') %% Control: non-target
    plot(epoav.t,epoav.x(:, i, 6),':b') %% Control: target
    title(chan{i}); axis(ax1);
    if i == 1
        legend(epoav.className);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  ============  3) Scalp Plot  =================================== 
clearvars -except path
cd(path)
load('epoav_includedSubjects.mat', 'epoav')
load('p3_range.mat')

% select interval
epoav  =  proc_selectIval(epoav, p3range);

% average of noMove and classic conditions
epoav1 = epoav;
epoav = proc_selectClasses(epoav, {'NoGo classic'  'NoGo syncMove'  'Go classic'  'Go syncMove'});
epoav.x = []; epoav.className = {'NoGo avNoMoveControl'  'NoGo syncMove'  'Go avNoMoveControl'  'Go syncMove'};
epoav.x(:, :, [2 4]) = epoav1.x(:, :, [2 5]);
epoav.x(:, :, 1) = mean(epoav1.x(:, :, [1 3]), 3); % noGo average between classic and noMove
epoav.x(:, :, 3) = mean(epoav1.x(:, :, [4 6]), 3); % go average between classic and noMove

% subtracting average - synchMove
epoav1 = epoav;
epoav = proc_selectClasses(epoav, {'NoGo classic' 'Go classic' });
epoav.x = []; epoav.className = {'NoGo difference'  'Go difference'};
epoav.x(:, :, 1) = epoav1.x(:, :, 1) - epoav1.x(:, :, 2);
epoav.x(:, :, 2) = epoav1.x(:, :, 3) - epoav1.x(:, :, 4);

% get electrode postitions
mnt = mnt_setElectrodePositions(epoav.clab);

% plot
figure;
for class = 1:2
    subplot(1, 2, class)
    plot_scalp(mnt, mean(epoav.x(:, :, class),1), 'CLim', [0 1.5]); title(epoav.className(class));
end
