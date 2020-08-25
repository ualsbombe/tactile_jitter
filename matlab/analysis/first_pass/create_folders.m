   
%% MEG
path =['/home/lau/projects/cerebellar_clock/scratch/tactile_jitter' ...
       '/MEG'];
subjects = {
    '0001' '0002' '0004' '0005' '0006' '0007' ...
    '0008' '0009' '0010' '0011' '0012' '0013' ...
    '0014' '0015' '0016' '0017' '0018' '0019' ...
    '0020' '0021' '0022' '0023' '0024' '0025' ...
    '0026' '0027' '0028' '0029' '0030' '0031'
    };
dates = {
    '20200124' '20200121' '20200121' '20200128' '20200120' '20200120' ...
    '20200120' '20200121' '20200122' '20200122' '20200122' '20200124' ...
    '20200127' '20200128' '20200128' '20200129' '20200129' '20200129' ...
    '20200131' '20200131' '20200131' '20200203' '20200203' '20200203' ...
    '20200204' '20200204' '20200204' '20200205' '20200205' '20200205'
         };
     
n_subjects = length(subjects);

for subject_index = 1:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    full_date = [date '_000000'];
    
    cd(path)
    mkdir(subject)
    cd(subject)
    mkdir(full_date)
    cd(path)
    
end

%% MR
path =['/home/lau/projects/cerebellar_clock/scratch/tactile_jitter' ...
       '/MEG'];
subjects = {
    '0001' '0002' '0004' ...
    '0005'
    };
dates = {
    '20180515_132117' '20191015_104445' '20191015_112257' ...
    '20191015_121553' 
        };
    
n_subjects = length(subjects);

for subject_index = 1:n_subjects
    
    subject = subjects{subject_index};
    date = dates{subject_index};
    
    cd(path)
    mkdir(subject)
    cd(subject)
    mkdir(date)
    cd(path)
    
end
     