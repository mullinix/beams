% plot_fra_curves
%%%%%%%%%%%%%%%%%%%%
% after looping over and solving the forced response, recall the data
% and generate plots of the forced response curves.

load matrices;
% list the files of the proper label.
% returns a long string of space delimited filenames.
% split the string by spaces, store in cell array.
phyles = strsplit(ls([pwd,'/data/forced_response_sweep/forced_response_data_real-*.mat']),'[\n]+');
% ignore "empty" strings from the split.
phyles=phyles(~strcmp(phyles,''));
% instantiate matrix of parameter values
params = zeros(length(phyles),2);
% initialize Y for plotting
Y=cell(length(phyles),1);
% for each file...
for idx=1:length(phyles)
    % load the data
    load(phyles{idx}); 
    % here we grab the value of gamma from the file name
    % separate name by the dash
    val_arry=strsplit(phyles{idx},'[-]+');
    % the third entry is the value of k
    kay = val_arry{3};
    % the end has the value of gamma and the extension
    gam=val_arry{end};
    % remove the extension
    gam=gam(1:end-4);
    %convert to a number
    gam=str2double(gam);
    kay=str2double(kay);
    params(idx,:)=[gam,kay];
    Y{idx}.Umax = max(abs(displacement.U));
    Y{idx}.Vmax = max(abs(displacement.V));
    Y{idx}.Wmax = max(abs(displacement.W));
end
[params,idx]=sortrows(params,[1,2]);
Y=Y(idx);
kk=7;
gam_idx=(params(:,2)>kk-0.1 & params(:,2)<kk+0.1);
gam=params(gam_idx,1);
Structs=cell2mat(Y(gam_idx));
Wplots=[Structs(:).Vmax];
plot(gam,Wplots(3:3:end));
% set(0,'defaulttextinterpreter','latex');
% figure(1);clf;
% hold on;
% % plot harmonic lines
% for i=1:10
% h1=plot(g./props.T*60,i*g./props.T*60,'--');
% end
% % plot omega vs gamma
% for idx=1:48
% h2=plot(g./props.T*60,Y(:,idx)./props.T*60);
% end
% xlim([0,max(g./props.T*60)]);ylim([0,3e5]);
% legendtext = {'$k_n=n\Omega_0$';'$\omega$ vs $\Omega_0$'};
% xlabel('$\Omega_0$','FontSize',15);
% ylabel('$\omega$','FontSize',15);
% % title({'Eigenfrequency $\omega$ vs'; 'hub rotation speed $\Omega_0$'},'FontSize',15)
% legend([h1,h2],legendtext,'interpreter','latex','FontSize',15,'location','northwest');