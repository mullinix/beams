#!/usr/bin/perl
if(@ARGV<1){
	$results = `cat git_results.txt`
}elsif(@ARGV>1) {
	die("Usage: # git_rm.pl <file in>"); 
}else {
	$results = `cat $ARGV[0]`;
}
@rm_arry = split /\n/,$results;
foreach (@rm_arry) {
	if( /deleted:/){
		@line_arry = split / +/,$_;
		$line_arry_size = @line_arry;
		print "$line_arry_size\n";
		if($line_arry_size>2) {
			$word_join_size = $line_arry_size-2;
			$file = join ' ',@line_arry[$word_join_size..$line_arry_size-1];
		}else{
			$file = $line_arry[@line_arry-1];
		}
		`git rm "$file"`;
		die("git failed to remove $file!\n$?\n") if($?);
	}
}

	
