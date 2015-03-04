%% get displacement modes
dof_shift = 5;
modes_size = num_beams*(num_nodes*node_dofs-num_bcs);
build_modes = zeros(floor(modes_size/5),1);
for i=1:num_beams
    length((i-1)*modes_size/num_beams/5+1:(i)*modes_size/num_beams/5)
    length((i-1)*modes_size/num_beams+dof_shift:5:i*modes_size/num_beams)
    build_modes((i-1)*modes_size/num_beams/5+1:(i)*modes_size/num_beams/5) = ...
        (i-1)*modes_size/num_beams+dof_shift:5:i*modes_size/num_beams;
end
modes1 = (shape(build_modes,:));% bending modes are in the last third
% modes1 = shape([dof_shift:5:end/2,end/2+dof_shift:5:end],:);% bending modes are in the last third
%% show that alternating cols are identical
% max(max(abs(modes1(:,2:2:end)-modes1(:,1:2:end-1))))
%% extract interesting non-duplicate modes (first 5 modes)
offset=1;
displ_modes = (modes1(:,(1:18)+offset-1));
displ_modes = real(displ_modes);
%% normalize to max of 1 (divide all elts by inf norm)
% for i=1:size(displ_modes,2)
%     displ_modes(:,i) = displ_modes(:,i)./max(abs(displ_modes(:,i)));
% end
%% or, normalize to total value of 1 (divide all elts by 2-norm)
% modes = normc(modes);
%% show what the shapes are
plot(displ_modes);
% %% get bending modes
% modes2 = shape(3:5:end/2,:);
% %% show that alternating cols are identical
% % sum(abs(modes2(:,2:2:end)-modes2(:,1:2:end)))
% bend_modes = modes2(:,1:2:10);
% %% normalize to max of 1 (divide all elts by inf norm)
% for i=1:size(bend_modes,2)
%     bend_modes(:,i) = bend_modes(:,i)./max(abs(bend_modes(:,i)));
% end
% figure();
% plot(bend_modes);