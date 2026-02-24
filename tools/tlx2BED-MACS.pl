#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use IO::File;
use Text::CSV;

my $tlxfile = shift;
my $bedfile = shift;
my $extsize = shift;

$extsize = 500 unless defined $extsize;

system("perl -pi -e 's/\r/\n/g' $tlxfile");

my $tlxfh = IO::File->new("<$tlxfile") or croak "Error: could not read file $tlxfile";
my $bedfh = IO::File->new(">$bedfile") or croak "Error: could not write file $bedfile";

my $csv = Text::CSV->new({sep_char=>"\t"});

my $header = $csv->getline($tlxfh);
$csv->column_names(map {lc} @$header);

my $i = 1;
while (my $read = $csv->getline_hr($tlxfh)) {
  my $rname = $read->{rname};
  my $rstart = $read->{strand} == 1 ? $read->{junction} - $extsize : $read->{junction} - 1;
  my $rend = $read->{strand} == 1 ? $read->{junction} : $read->{junction} + $extsize - 1;
  my $strand = $read->{strand} == 1 ? "-" : "+";

  $bedfh->print(join("\t",$rname,$rstart,$rend,"TLX_".$i,0,$strand)."\n");
  $i++;
}