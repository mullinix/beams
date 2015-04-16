%% forced response loop
% list the files of the proper label.
% returns a long string of space delimited filenames.
% split the string by spaces, store in cell array.
phyles = strsplit(ls([pwd,'/data/omega_evals_sweep/datastruct-*.mat']),'[\n]+');
% ignore "empty" strings from the split.
phyles=phyles(~strcmp(phyles,''));
% instantiate vecror of gamma values
g = zeros(length(phyles),1);
% initialize Y for plotting
% Y=zeros(length(phyles),1);
% for each file...
for idx=1:length(phyles)
    % load the data
    load(phyles{idx}); 
    % here we grab the value of gamma from the file name
    % separate name by the dash
    gam=strsplit(phyles{idx},'[-]+');
    % the end has the value of gamma and the extension
    gam=gam{end};
    % remove the extension
    gam=gam(1:end-4);
    %convert to a number
    gam=str2double(gam);
    g(idx)=gam;
    freq = data2save.freqs(1);
    forced_response;
    
end