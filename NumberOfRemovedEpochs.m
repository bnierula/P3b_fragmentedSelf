%% calculate number of removed channels and epochs to report it in the methods part of the paper

rootPath = 'E:\DATA_EventLab\Projects\Noself3Pb\Analysis\EEG\';

nChan = zeros(35,1);
chanName = cell(35,1);
nEpo = NaN(35,1);

for experiment=1:2  

    if experiment==1
        myPath=[rootPath 'ExpJune\'];
        nSubjects=[1:2 4:25];
    elseif experiment==2
        myPath=[rootPath 'ExpOctober\'];
        nSubjects=26:35;
    end
    cd(myPath)

    for subject=nSubjects
        if subject<10
            cd([myPath 'S0' num2str(subject)])
        else
            cd([myPath 'S' num2str(subject)])
        end
        
        load('artifacts.mat', 'art')
        
%         if ~isempty(art.clabNot)
%             if isempty(art.clabNot)
%                 nChan(subject,1) = 0;
%             else
%                 nChan(subject,1) = length(art.clabNot);
%             end
%         end
%         chanName(subject,1) = {art.clabNot};
        nEpo(subject,1) = length(art.epoNot);
    end
end

% statistics on excluded epochs
meanEpo = mean(nEpo([1 4 5 7 9:21 23 25 27:29 31:35], 1));
stdEpo = std(nEpo([1 4 5 7 9:21 23 25 27:29 31:35], 1));

% statistics on excluded channels
% sumChan = sum(nChan);