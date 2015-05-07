% generate avis of forced response displacement data.

% list files to make movies from
phyles=strsplit(ls([pwd,'/data/forced_response_sweep/bu2/forced_response_data_real*']),'[\n]+');
% ignore "empty" strings from the split.
phyles=phyles(~strcmp(phyles,''));
for idx=1:length(phyles)
    % load the data
    load(phyles{idx}); 
    freq=strsplit(phyles{idx},'[-]+');
    % the second entry has the value of freq
    freq=freq{3};
%     load('./data/forced_response_sweep/forced_response_data_real-freq-5.9-gamma-5.0.mat')
%     freq='5.9';
    moov(1:length(displacement.U)) = struct('cdata',[],'colormap',[]);

    circ_x = cos(0:0.01:2*pi)*props.a;
    circ_y = sin(0:0.01:2*pi)*props.a;
    circ_z = 0.*(0:0.01:2*pi);
    scale=1e2;
    lsz=3;
    mks=25;

    x_len = linspace(0,1,length(displacement.U(1,:))/3);
    y_len = linspace(0,1,length(displacement.V(1,:))/3);
    beam_len = length(x_len);

    for i=1:length(displacement.U)
    figure(1); clf;
    plot3(circ_x,circ_y,circ_z,'LineWidth',lsz);
    hold on;

    xacks = (props.L+props.a+max(max(abs(displacement.U.*scale))))*1.05;
    yacks = (props.L+props.a+max(max(abs(displacement.V.*scale))))*1.05;
    macks2d=max([xacks,yacks])*.75;
    zacks = (max(max(abs(displacement.W.*scale))))*1.05+1e-16;

    for j=1:num_beams
        plt_pts = (j-1)*beam_len+1:j*beam_len;
        x=cos(angles(j)).*(props.L.*x_len+props.a);
        y=sin(angles(j)).*(props.L.*y_len+props.a);
        plot3(x+displacement.U(i,plt_pts).*scale,...
            y+displacement.V(i,plt_pts).*scale,...
            displacement.W(i,plt_pts).*scale,'LineWidth',lsz);
        plot3(x(end)+displacement.U(i,plt_pts(end)).*scale,...
            y(end)+displacement.V(i,plt_pts(end)).*scale,...
            displacement.W(i,plt_pts(end)).*scale,'r.','MarkerSize',mks);
        xlim([-macks2d,macks2d]);
        ylim([-macks2d,macks2d]);
        zlim([-zacks,zacks]);
    end
    moov(i)=getframe(); %#ok<SAGROW>

    end
    movie2avi(moov, sprintf('%s/data/forced_response_sweep/freq-%s.avi',pwd,freq),...
        'compression', 'None');
end
