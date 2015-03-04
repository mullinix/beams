build_coupling_stiffness;
evec_perm;
z=zeros(size(shape(:,mode)));
Y0 = [shape(:,mode);
                  z; ];
numpts = 1000;
num_cycles = 2;
t_final = (2*pi*num_cycles)/(freqs(mode));
t = linspace(0,t_final,numpts);
X2inv = inv(X2);
X1 = X1*1i;
%% ode solution
[t,Y] = ode23t(@(t,Y)beam_ode(t,Y,X0,X1,X2inv),t,Y0);

beampts = length(Y(1,1:end/2));
beamx = linspace(0,1,beampts/num_beams);
beams_x = ones(beampts,1);
beams_y = ones(beampts,1);
for i=1:length(angles)
    beams_x((i-1)*floor(beampts/num_beams)+1:i*floor(beampts/num_beams)) = cos(angles(i)).*(props.L.*beamx+props.a);
    beams_y((i-1)*floor(beampts/num_beams)+1:i*floor(beampts/num_beams)) = sin(angles(i)).*(props.L.*beamx+props.a);
end
circ_x = cos(0:0.01:2*pi)*props.a;
circ_y = sin(0:0.01:2*pi)*props.a;
circ_z = 0.*(0:0.01:2*pi);
xmax = max(abs(beams_x))*1.25;
ymax = max(abs(beams_y))*1.25;
% xmax = max([xmax,ymax,abs(circ_x)])*1.25;
% ymax = xmax;%max([xmax,ymax,abs(circ_x)])*1.25;
start_pt = 5;

high = max(max(abs(Y(:,start_pt:5:end))));
high = high*1.25;

%% animate
for i=1:numpts
    figure(1);clf;
    hold off;
    plot3(circ_x,circ_y,circ_z);
    hold on;
    for j=1:num_beams
        beam_pts = beampts*(j-1)/num_beams+start_pt:5:(j)*beampts/num_beams;
        plot3(beams_x(beam_pts),beams_y(beam_pts),real(Y(i,beam_pts)));
    end

    xlim([-xmax,xmax]);
    ylim([-ymax,ymax]);
    zlim([-high,high]);
    drawnow;
end