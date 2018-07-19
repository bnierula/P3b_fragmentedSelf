%% prepare variables to copy into excel
myPath = 'E:\DATA_Eventlab\Projects\Noself3Pb\Analysis\EEG\';
addpath(genpath([myPath 'Scripts']))
cd(myPath)

load('p3.mat', 'p3b'); % 1=NoGo classic; 2=NoGo syncMove; 3=NoGo noMove; 
                       % 4=Go classic; 5=Go syncMove; 6=Go noMove
p3 = reorder4Statistics(p3b);     
clear p3b

load('n1.mat', 'po3'); 
n1po3 = reorder4Statistics(po3); 
clear po3

load('n1.mat', 'po4'); 
n1po4 = reorder4Statistics(po4); 
clear po4

load('p1.mat', 'po3'); 
p1po3 = reorder4Statistics(po3); 
clear po3

load('p1.mat', 'po4'); 
p1po4 = reorder4Statistics(po4); 
clear po4

