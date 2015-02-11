%% get displacement modes
% modes1 = shape(2*end/3+3:5:end,:);% bending modes are in the last third
%% show that alternating cols are identical
max(max(abs(modes1(:,2:2:end)-modes1(:,1:2:end-1))))
%% extract interesting non-duplicate modes (first 5 modes)
displ_modes = real(modes1(:,end:-1:end-9));
%% normalize to max of 1 (divide all elts by inf norm)
for i=1:size(displ_modes,2)
    displ_modes(:,i) = displ_modes(:,i)./max(abs(displ_modes(:,i)));
end
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