function splitstr = strsplit(str_in,split_in)

[~, splitstr] = regexp(str_in, split_in, 'match', 'split');

end