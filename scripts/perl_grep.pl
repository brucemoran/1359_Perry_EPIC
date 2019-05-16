#!/usr/bin/perl -w

##https://stackoverflow.com/questions/42239179/fastest-way-to-find-lines-of-a-file-from-another-larger-file-in-bash

##add parseheaderkey to select columns by keyword
use strict;
if (scalar(@ARGV) != 3) {
  printf STDERR "Usage: fgrep.pl smallfile bigfile parseheaderkey\n";
  exit(2);
}

my ($small_file, $big_file) = ($ARGV[0], $ARGV[1]);
my ($small_fp, $big_fp, %small_hash, @field, @cols);

open($small_fp, "<", $small_file) || die "Can't open $small_file: " . $!;
open($big_fp, "<", $big_file)     || die "Can't open $big_file: "   . $!;

# store contents of small file in a hash
while (<$small_fp>) {
  chomp;
  $small_hash{$_} = 1;
}
close($small_fp);

# loop through big file and find matches
my $a=0;
while (<$big_fp>) {

  ##split header, find cols to output
  if($a==0){
    my @index=split(/\t/, $_);
    push(@cols,0);
    print "CpGs";
    for(my $i=0;$i<@index;$i++){
      if($index[$i]=~m/$ARGV[2]/){
        push(@cols,$i);
        print "," . $index[$i];
      }
    }
    print "\n";
    $a++;
    next;
  }

  @field=split(/\t/, $_);
  if(defined($field[0]) && exists($small_hash{$field[0]})) {
    my @out;
    foreach my $o (@cols){
      push(@out,$field[$o]);
    }
    print join(",", @out) . "\n";
  }
}

close($big_fp);
exit(0);
