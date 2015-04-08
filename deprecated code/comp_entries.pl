#!/usr/bin/perl
$N = @ARGV;
if($N!=1) {
	print "Usage: $ ./compare_entries.pl <input_file>\n";
	exit(0);
} else {
	$fname = $ARGV[0];
}
open(SEARCH,$fname);
%copy_strs = ();
%copy_entries = ();
$j=0;
while(<SEARCH>) {
	chomp();
	$line = $_;
	$j++;
	@line_arry = split(/ /,$line);
	$first_val = $line_arry[1];
	$first_val =~ s/^\s+|\s+$//g;
	if( $first_val ne "0" ) {
		$entry = $line_arry[0];
		$idx = 1;
		$arry_size = @line_arry;
		$line_str = join(" ",@line_arry[1..$#line_arry]);
		$uniq=1;
		$line_str =~ s/^\s+|\s+$//g;
		foreach $key (sort(keys %copy_strs)) {
			$value = $copy_strs{$key};
			if($line_str eq $value) {
				push(@{$copy_entries{$key}}, $entry);
				$uniq=0;
				last;
			}
			if($line_str eq "-".$value) {
				push(@{$copy_entries{$key}}, "-".$entry);
				$uniq=0;
				last;
			}
		}
		if($uniq==1) {
			$copy_strs{$entry} = $line_str;
		}
	}
}
close(SEARCH);
$i=1;
$total_unique = keys %copy_strs;
print "totes: $total_unique\n";
foreach $key (sort(keys %copy_strs)) {
	$value = $copy_strs{$key};
	print "$i $key $value\n\t";
	foreach $copy (@{$copy_entries{$key}}) {
		print " $copy ";
	}
	print "\n";
	$i++;
}
#eof

