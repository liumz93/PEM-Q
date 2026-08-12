#!/usr/bin/env perl


use strict;
use warnings;
use Getopt::Long;
use Carp;
use IO::Handle;
use IO::File;
use Text::CSV;
use File::Basename;
use File::Which;
use List::Util qw(min max);
use Interpolation 'arg:@->$' => \&argument;
use Time::HiRes qw(gettimeofday tv_interval);

use Cwd qw(abs_path);
use FindBin;
use lib abs_path("$FindBin::Bin/../lib");

require "PerlSub.pl";
require "TranslocHelper.pl";


# Flush output after every write
select( (select(STDOUT), $| = 1 )[0] );

##
## This program
## 
## run with "--help" for usage information
##
## Robin Meyers

# Forward declarations
sub parse_command_line;
sub read_in_translocations;
sub align_to_primer ($);
sub align_to_adapter ($);
sub print_html_header ($);
sub print_html_footer ($);
sub write_marked_read ($$);


# Global flags and arguments, 
# Set by command line arguments
my $tlxfile;
my $htmlfile;
my $primer = "";
my $adapter = "";


# Global variables
my $htmlfh;
my $tlxfh;

#
# Start of Program
#

parse_command_line;

my $t0 = [gettimeofday];

$htmlfh = IO::File->new(">$htmlfile") or croak "Error: cannot write to html file $htmlfile";
print_html_header($htmlfh);


$tlxfh = IO::File->new("<$tlxfile");
my $csv = Text::CSV->new({sep_char => "\t"});
my $header = $csv->getline($tlxfh);
$csv->column_names(@$header);


while (my $tl = $csv->getline_hr($tlxfh)) {

  align_to_primer($tl) if $primer =~ /\S/;
  align_to_adapter($tl) if $adapter =~ /\S/;

  write_marked_read($htmlfh, $tl);

}

$tlxfh->close;

print_html_footer($htmlfh);


$htmlfh->close;



my $t1 = tv_interval($t0);

printf("\nFinished writing HTML reads in %.2f seconds.\n", $t1);



#
# End of program
#

sub write_marked_read ($$) {
  my $fh = shift;
  my $tl = shift;

  my @seq = split("",$tl->{Sequence});
  my %marks;
  my $marktypes = {alignment => {start => "Qstart", end => "Qend"},
                   breaksite => {start => "B_Qstart", end => "B_Qend"},
                   primer => {start => "PrimQstart", end => "PrimQend"},
                   adapter => {start => "AdptQstart", end => "AdptQend"}};

  foreach my $marktype (sort keys %$marktypes) {

    if (exists $tl->{$marktypes->{$marktype}->{start}}) {

      my $start = $tl->{$marktypes->{$marktype}->{start}} - 1;
      my $end = $tl->{$marktypes->{$marktype}->{end}} - 1;

      if (exists $marks{$start}) {
        push(@{$marks{$start}->{start}},$marktype);
      } else {
        $marks{$start}->{start} = [$marktype];
      }

      if (exists $marks{$end}) {
        push(@{$marks{$end}->{end}},$marktype);
      } else {
        $marks{$end}->{end} = [$marktype];
      }
    }
  }

  my %open_marks;
  my $marked_read = "<div class=\"sequence\">";

  foreach my $i (0..$#seq) {

    if (exists $marks{$i}->{start}) {

      my @marks = @{$marks{$i}->{start}};

      $marked_read .= "</span>" if keys %open_marks > 0;
      $marked_read .= "<span class=\"" . join(" ",@marks) . 
                    (keys %open_marks > 0 ? " ".join(" ",keys %open_marks) : "") . "\">";

      foreach my $mark (@marks) {
        $open_marks{$mark} = 1;
      }
    }

    $marked_read .= $seq[$i];

    if (exists $marks{$i}->{end}) {
      my @marks = @{$marks{$i}->{end}};

      $marked_read .= "</span>";

      foreach my $mark (@marks) {
        delete $open_marks{$mark};
      }

      $marked_read .= "<span class=\"" . join(" ",keys %open_marks) . "\">" if keys %open_marks > 0;

    }

  }

  $marked_read .= "</span>" if keys %open_marks > 0;

  $marked_read .= "</div>";

  my $qname = $tl->{Qname};
  my $rname = $tl->{Rname};
  my $rstart = $tl->{Prey_start};
  my $rend = $tl->{Prey_end};
  my $strand = $tl->{Strand};

  my $marked_id = "<div class=\"seqname\">>$qname $rname:$rstart-$rend:$strand</div>";
  $fh->print("$marked_id\n");
  $fh->print("$marked_read\n");


}


sub align_to_primer ($) {
	
  my $tl = shift;
  
  my $strand = $tl->{Bait_strand};
  
  # my $primer = reverse($primer) if ($strand eq "-");
  
  my $qstart = index($tl->{Sequence}, $primer);

  unless ($qstart == -1) {
    $tl->{PrimQstart} = $qstart + 1;
    $tl->{PrimQend} = $qstart + length($primer);
  }

}

sub align_to_adapter ($) {
	
  my $tl = shift;
  
  my $strand = $tl->{Bait_strand};
  # my $adapter = reverse($adapter) if ($strand eq "-");

  my $qstart = index($tl->{Sequence}, $adapter);

  unless ($qstart == -1) {
    $tl->{AdptQstart} = $qstart + 1;
    $tl->{AdptQend} = $qstart + length($adapter);
  }
}



sub parse_command_line {
  my $help;

  usage() if (scalar @ARGV == 0);

  my $result = GetOptions ( "primer=s" => \$primer,
                            "adapter=s" => \$adapter,
                            "help" => \$help
                          ) ;
  
  usage() if ($help);

  $tlxfile = shift(@ARGV);
  $htmlfile = shift(@ARGV);

  #Check options

  croak "Error: cannot read tlx file $tlxfile" unless -r $tlxfile;
  
  exit unless $result;
}

sub print_html_header ($) {
  my $fh = shift;
  print $fh <<EOF
<!DOCTYPE html>
<html>
<head>
  <style type="text/css">
      #sequences {
        word-wrap:break-word;
      }
      .sequence {
        letter-spacing:-1px;
        font-family:courier;
      }
      .primer {
        color:red;
      }
      .adapter {
        color:blue;
      }
      .breaksite {
        text-decoration:underline;
      }
      .alignment {
        background-color:yellow;
      }
  </style>
</head>

<body>
<div id="sequences">
EOF
}

sub print_html_footer ($) {
  my $fh = shift;
  print $fh <<EOF
</div>
</body>
</html>
EOF
}


sub usage()
{
print<<EOF;
Title, by Robin Meyers, ddmonthyyyy

This program .


Usage: $0 arg1 arg2 arg3 ...
        [--option VAL] [--flag] [--help]

Arguments (defaults in parentheses):

$arg{"tlxfile"," "}
$arg{"htmlfile"," "}
$arg{"--primer"," ",$primer}
$arg{"--adapter"," ",$adapter}
$arg{"--help","This helpful help screen."}


EOF

exit 1;
}
