% find_modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% will plot several shapes for a mode to visualize the system.

%% get displacement modes
dof_shift = 5;
chunk_size = 2;
plotEvals = SVevals;
plotShape = SVshape;
elts = grab_elts(chunk_size,length(plotEvals)/num_beams/2,1,dof_shift);
build_modes = repmat(elts,1,3);
modes1 = (plotShape(build_modes,:));% bending modes are in the last third
%% show that alternating cols are identical
% max(max(abs(modes1(:,2:2:end)-modes1(:,1:2:end-1))))
%% extract interesting non-duplicate modes (first 5 modes)
offset=1;
displ_modes = (modes1(:,(1:10)+offset-1));
displ_modes = real(displ_modes);
%% normalize to max of 1 (divide all elts by inf norm)
for i=1:size(displ_modes,2)
    displ_modes(:,i) = displ_modes(:,i)./max(abs(displ_modes(:,i)));
end
%% or, normalize to total value of 1 (divide all elts by 2-norm)
% modes = normc(modes);
%% show what the shapes are
plot(displ_modes);
