function cnt1 = interpolateChannels(cnt, chanIdx)

%% save data
path1 = pwd;
path = strrep(path1,'/','//');
X = cnt.x;
save temp X
clear X


%% Add path for eeglab toolbox
% addpath(eeglabPath);

%% start eeglab
eeglab

%% import data from matlab file
EEG.etc.eeglabvers = '14.1.1'; % this tracks which version of EEGLAB is being used, you may ignore it
data = [ path '//temp.mat'];
setname = cnt.title;
srate = cnt.fs;
subject = cnt.title(11:12);
pnts = 0;
condition = cnt.title(1:9);
xmin = 0;
session = 1;
idcs = strfind(path1,filesep);
chanlocs = [path1(1:idcs(9)) 'Scripts/channelLocations.sfp'];
chanlocs = strrep(chanlocs,'/','//');
EEG = pop_importdata('dataformat','matlab','nbchan',0,'data',data,...
    'setname',setname,'srate',srate,'subject',subject,'pnts',pnts,...
    'condition',condition,'xmin',xmin,'session',session,'chanlocs',chanlocs);
EEG = eeg_checkset( EEG );

%% interpolate channels
EEG = pop_interp(EEG, chanIdx, 'spherical');

%% transform back to BBCI toolbox format
cnt1 = cnt; cnt1.x = [];
cnt1.x = double(EEG.data');

delete([pwd '/temp.mat'])

%% plot
for i=1:length(chanIdx)
    figure; 
    plot(cnt.x(:,chanIdx(i)),'b'); hold on;
    plot(cnt1.x(:,chanIdx(i)),'r');
    pause(5)
    close
end


%% stop eeglab by code
close 