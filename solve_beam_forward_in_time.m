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
X1 = -X1*1i;
%% ode solution
[t,Y] = ode23t(@(t,Y)beam_ode(t,Y,X0,X1,X2inv),t,Y0);

% beampts = length(Y(1,1:end/2));
% beamx = linspace(0,1,beampts/num_beams);
% beams_x = ones(beampts,1);
% beams_y = ones(beampts,1);
% for i=1:length(angles)
%     beams_x((i-1)*floor(beampts/num_beams)+1:i*floor(beampts/num_beams)) = cos(angles(i)).*(props.L.*beamx+props.a);
%     beams_y((i-1)*floor(beampts/num_beams)+1:i*floor(beampts/num_beams)) = sin(angles(i)).*(props.L.*beamx+props.a);
% end
% xmax = max(abs(beams_x))*1.25;
% ymax = max(abs(beams_y))*1.25;
% xmax = max([xmax,ymax,abs(circ_x)])*1.25;
% ymax = xmax;%max([xmax,ymax,abs(circ_x)])*1.25;

circ_x = cos(0:0.01:2*pi)*props.a;
circ_y = sin(0:0.01:2*pi)*props.a;
circ_z = 0.*(0:0.01:2*pi);

start_pt = 5;

high = max(max(abs(Y(:,start_pt:5:end))));
high = high*1.25;

spts = grab_elts(5,length(evals)/num_beams/2,1,2);
spts = repmat(spts,1,3);
s_max = max(max(abs(Y(:,spts))));
v_pts = grab_elts(5,length(evals)/num_beams/2,1,3);
v_pts = repmat(v_pts,1,3);
v_max = max(max(abs(Y(:,v_pts))));
w_pts = grab_elts(5,length(evals)/num_beams/2,1,5);
w_pts = repmat(w_pts,1,3);
w_max = max(max(abs(Y(:,w_pts))));
v_full = grab_elts(5,length(evals)/num_beams/2,2,3);
v_full = repmat(v_full,1,3);
w_full = grab_elts(5,length(evals)/num_beams/2,2,5);
w_full(1) = (1==0);
w_full = repmat(w_full,1,3);
x=linspace(0,props.L,num_elts+1);
x=x(2:end);

% plotting axis sizes
xmax = (max(abs((props.L+props.a).*cos(angles)))+s_max)*1.1;
ymax = (max(abs((props.L+props.a).*sin(angles)))+v_max)*1.1;
zmax = w_max*1.1;
if(xmax<1e-3)
    xmax=1e-3;
end
if(ymax<1e-3)
    ymax=1e-3;
end
if(zmax<1e-3)
    zmax=1e-3;
end

%% animate
for i=1:numpts
    Si = Y(i,spts);  %Si = Si*5;
    Vi = Y(i,v_pts); %Vi = Vi*5;
    Wi = Y(i,w_pts); %Wi = Wi*5;
    Vfull = Y(i,v_full);
    Wfull = Y(i,w_full);
    len=length(Vi);
    len2=length(Vfull);
    figure(1);clf;
    hold off;
    plot3(circ_x,circ_y,circ_z);
    hold on;
    for j=1:num_beams
        plot_pts = (len/num_beams*(j-1)+1):len/num_beams*j;
        plot_pts2 = (len2/num_beams*(j-1)+1):len2/num_beams*j;
        H_vw = hvw(x,Vfull(plot_pts2),Wfull(plot_pts2));
        Ui = Si(plot_pts)-H_vw;
        u=real(Ui);           u = u.*5;
        v=real(Vi(plot_pts)); v = v.*5;
        w=real(Wi(plot_pts)); w = w.*5;
        x_len = linspace(0,1,length(u));
        y_len = linspace(0,1,length(v));
        x=cos(angles(j)).*(props.L.*x_len+props.a);
        y=sin(angles(j)).*(props.L.*y_len+props.a);
        plot3(x+u,y+v,w);
    end
    xlim([-xmax,xmax]);
    ylim([-ymax,ymax]);
    zlim([-zmax,zmax]);
%     view(2);
    drawnow;
end