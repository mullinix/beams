#!/usr/bin/perl
if(@ARGV<1){
	$results = `cat git_results.txt`
}elsif(@ARGV>1) {
	die("Usage: # git_rm.pl <file in>"); 
}else {
	$results = `cat $ARGV[0]`;
}
#print $results;
@rm_arry = split /\n/,$results;
print "$rm_arry[0]\n";
foreach (@rm_arry) {
	if( /deleted:/){
		@line_arry = split / /,$_;
		$file = $line_arry[@line_arry-1];
		`git rm $file`;
		die("git failed to remove $file!\n$?\n") if($?);
#		print "$file\n";
	}
}

	
