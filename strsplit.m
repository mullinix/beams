function splitstr = strsplit(str_in,split_in)
% function splitstr = strsplit(str_in,split_in)
% Returns an array of strings split by "split_in"

[~, splitstr] = regexp(str_in, split_in, 'match', 'split');

end