% plot_gamma_vs_omega
%%%%%%%%%%%%%%%%%%%%%%%
% generates plots of intersections of n*Omega and the eigenvalue curves

load matrices;
% list the files of the proper label.
% returns a long string of space delimited filenames.
% split the string by spaces, store in cell array.
phyles = strsplit(ls([pwd,'/data/omega_evals_sweep/datastruct-*.mat']),'[\n]+');
% ignore "empty" strings from the split.
phyles=phyles(~strcmp(phyles,''));
% instantiate vecror of gamma values
g = zeros(length(phyles),1);
% initialize Y for plotting
Y=zeros(length(phyles),48);
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
    Y(idx,:) = data2save.omega;
end
[g,idx]=sort(g);
Y = Y(idx,:);
set(0,'defaulttextinterpreter','latex');
figure(1);clf;
hold on;
% plot harmonic lines
for i=1:10
h1=plot(g./props.T*60,i*g./props.T*60,'k--');
% h1=plot(g,i*g,'k--');

end
% plot omega vs gamma
for idx=1:48
h2=plot(g./props.T*60,Y(:,idx)./props.T*60);
% h2=plot(g,Y(:,idx));
end
xlim([0,max(g./props.T*60)]);ylim([0,3e5]);
legendtext = {'$k_n=n\Omega_0$';'$\omega$ vs $\Omega_0$'};
xlabel('$\Omega_0$','FontSize',15);
ylabel('$\omega$','FontSize',15);
% title({'Eigenfrequency $\omega$ vs'; 'hub rotation speed $\Omega_0$'},'FontSize',15)
legend([h1,h2],legendtext,'interpreter','latex','FontSize',15,'location','northwest');