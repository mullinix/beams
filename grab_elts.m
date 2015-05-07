function elts = grab_elts(chunk_size, array_length, group_size, idx_start) 
%
%     function elts = grab_elts(chunk_size, array_length, group_size, idx_start)
%     A function to grab sub-elements. James Mullinix -- 3/11/15
%     
%     inputs:
%         chunk_size = (int) Size of repeating groups. (i.e. the matrix has
%             data in chunks of 5)
%         array_length = (int) Size of the vector (or matrix)
%         group_size = (int) Number of sequential elements to grab (i.e.
%             chunk size is 5, say, and we want groups of 3)
%         idx_start = (int) Index to start the group grab. (i.e. chunk size 
%             is 5, say, and we want groups of 3, starting at the second
%             element, set this to 2)
% 
%     output:
%         elts = (logical array) One if we grab the element, zero otherwise.
% 
%     Hint: If you need asymmetrical data grab, say in our case of 5 elts, we
%     need groups of 2, but we need to have element 1, skipping 2-4, grabbing
%     5 and 6, etc., then set idx_start to 5.
%     
%     Ref: http://www.mathworks.com/matlabcentral/answers/80246-i-need-to-remove-every-other-9-rows-of-a-matrix
%
modtester = mod((1:array_length)-1,chunk_size)+1;
idx=mod((idx_start:(idx_start+group_size-1))-1,chunk_size)+1;
elts=ones(size(modtester))==0;
for i=1:length(idx)
    elts=elts|modtester==idx(i);
end

end

