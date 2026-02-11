#!/usr/bin/perl

use strict;
use warnings;
use Carp;
use File::Basename;
use IO::Handle;
use IO::File;
use CGI;
use List::Util qw(sum);
use IPC::System::Simple qw(system capture);


# from Perl Cookbook http://docstore.mik.ua/orelly/perl/cookbook/ch08_09.htm
# usage: build_index(*DATA_HANDLE, *INDEX_HANDLE)
sub build_index {
	my $data_file  = shift;
	my $index_file = shift;
	my $offset     = 0;

	while (<$data_file>) {
		print $index_file pack("N", $offset);
		$offset = tell($data_file);
	}
}
# from Perl Cookbook http://docstore.mik.ua/orelly/perl/cookbook/ch08_09.htm
# usage: line_with_index(*DATA_HANDLE, *INDEX_HANDLE, $LINE_NUMBER)
# returns line or undef if LINE_NUMBER was out of range
sub line_with_index {
	my $data_file   = shift;
	my $index_file  = shift;
	my $line_number = shift;

	my $size;               # size of an index entry
	my $i_offset;           # offset into the index of the entry
	my $entry;              # index entry
	my $d_offset;           # offset into the data file

	$size = length(pack("N", 0));
	$i_offset = $size * ($line_number-1);
	seek($index_file, $i_offset, 0) or return;
	read($index_file, $entry, $size);
	$d_offset = unpack("N", $entry);
	seek($data_file, $d_offset, 0);
	return scalar(<$data_file>);
}

# from Perl Cookbook http://docstore.mik.ua/orelly/perl/cookbook/ch08_09.htm
# usage: line_with_index(*DATA_HANDLE, *INDEX_HANDLE, $LINE_NUMBER)
# returns line or undef if LINE_NUMBER was out of range
sub line_with_index_hr {
	my $data_file   = shift;
	my $csv_obj     = shift;
	my $index_file  = shift;
	my $line_number = shift;

	my $size;               # size of an index entry
	my $i_offset;           # offset into the index of the entry
	my $entry;              # index entry
	my $d_offset;           # offset into the data file

	$size = length(pack("N", 0));
	$i_offset = $size * ($line_number-1);
	seek($index_file, $i_offset, 0) or return;
	read($index_file, $entry, $size);
	$d_offset = unpack("N", $entry);
	seek($data_file, $d_offset, 0);
	return scalar($csv_obj->getline_hr($data_file));
}

sub mean {
	return sum(@_)/@_;
}

sub argument {
	my $var = shift;
	my $description = shift;
	my $default = shift;

	return sprintf("  \%\-16s - %s\n",
		(defined $default ? sprintf("%s (%s)",$var,$default) : sprintf("%s",$var)),
		$description);
}

sub read_fastq ($) {

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = $fh->getline();
	return () unless defined $seq_name;
	chomp($seq_name);
	croak "Error: unexpected format - missing '\@' symbol in sequence header $seq_name" unless $seq_name =~ s/^\@//;
	my $seq_bases = $fh->getline();
	croak "Error: bad input file, expecting line with sequences" unless defined $seq_bases;
	croak "Error: unexpected base in sequence $seq_bases" unless $seq_bases =~ /^[AGCTagctNn]*$/;
	chomp($seq_bases);
	my $seq_name2 = $fh->getline();
	croak "Error: bad input file, expecting line with sequence name2" unless defined $seq_name2;
	croak "Error: unexpected format - missing '\+' symbol in sequence name2 $seq_name2" unless $seq_name2 =~ /^\+.*$/;
	chomp($seq_name2);
	my $seq_quals = $fh->getline();
	croak "Error: bad input file, expecting line with quality scores" unless defined $seq_quals;
	chomp($seq_quals);
	croak "Error: sequence and quality must be same length" unless length($seq_bases) == length($seq_quals);

	return ($seq_name,$seq_bases,$seq_name2,$seq_quals);
}

sub write_fastq ($$$$$) {
	croak "Error: too few arguments to 'write_fastq' method" if scalar @_ < 5;

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = shift;
	my $seq_bases = shift;
	my $seq_name2 = shift;
	my $seq_quals = shift;

	$seq_name = "@" . $seq_name unless $seq_name =~ /^\@.*$/;
	croak "Error: sequence contains unexpected base" unless $seq_bases =~ /^[AGCTagctNn]+$/;
	$seq_name2 = "+" . $seq_name2 unless $seq_name2 =~ /^\+.*$/;
	croak "Error: sequence and quality must be same length" unless length($seq_bases) == length($seq_quals);

	print $fh $seq_name,"\n";
	print $fh $seq_bases,"\n";
	print $fh $seq_name2,"\n";
	print $fh $seq_quals,"\n";

}

sub read_fasta ($) {

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = $fh->getline();
	return () unless defined $seq_name;
	chomp($seq_name);
	croak "Error: unexpected format - missing '\>' symbol in sequence header $seq_name" unless $seq_name =~ s/^\>\s*(\S.*)$/$1/;
	my $seq_bases = "";
	while (my $line = $fh->getline()) {
		if ($line =~ /^>/) {
			seek($fh,-length($line),1);
			last;
		}
		chomp($line);
		croak "Error: unexpected base in sequence $line" unless $line =~ /^[AGCTagctNn\-\.]*$/;
		$seq_bases .= $line;
	}
	
	return ($seq_name,$seq_bases);
}

sub write_fasta ($$$) {
	croak "Error: too few arguments to 'write_fastq' method" if scalar @_ < 3;

	my $fh = shift;
	croak "Error: invalid file handle" unless $fh->opened;
	my $seq_name = shift;
	my $seq_bases = shift;
 

	$seq_name = ">" . $seq_name unless $seq_name =~ /^\>.*$/;
	croak "Error: sequence contains unexpected base" unless $seq_bases =~ /^[AGCTagctNn\-\.]+$/;

	print $fh $seq_name,"\n";
	print $fh $seq_bases,"\n";


}

sub blast_header {
	return qw(qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore);
}

sub quals_to_ascii (@) {
	my @quals = @_;
	return join("",map( chr($_ + 33), @quals ));
}

sub quals_from_ascii ($) {
	my $quals = shift;
	return map( ord($_) - 33, split("", $quals));
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }



# Return the Levenshtein distance (also called Edit distance) 
# between two strings
#
# The Levenshtein distance (LD) is a measure of similarity between two
# strings, denoted here by s1 and s2. The distance is the number of
# deletions, insertions or substitutions required to transform s1 into
# s2. The greater the distance, the more different the strings are.
#
# The algorithm employs a promimity matrix, which denotes the 
# distances between substrings of the two given strings. Read the 
# embedded comments for more info. If you want a deep understanding 
# of the algorithm, printthe matrix for some test strings 
# and study it
#
# The beauty of this system is that nothing is magical - the distance
# is intuitively understandable by humans
#
# The distance is named after the Russian scientist Vladimir
# Levenshtein, who devised the algorithm in 1965
#
sub levenshtein
{
    # $s1 and $s2 are the two strings
    # $len1 and $len2 are their respective lengths
    #
    my ($s1, $s2) = @_;
    my ($len1, $len2) = (length $s1, length $s2);

    # If one of the strings is empty, the distance is the length
    # of the other string
    #
    return $len2 if ($len1 == 0);
    return $len1 if ($len2 == 0);

    my %mat;

    # Init the distance matrix
    #
    # The first row to 0..$len1
    # The first column to 0..$len2
    # The rest to 0
    #
    # The first row and column are initialized so to denote distance
    # from the empty string
    #
    for (my $i = 0; $i <= $len1; ++$i)
    {
        for (my $j = 0; $j <= $len2; ++$j)
        {
            $mat{$i}{$j} = 0;
            $mat{0}{$j} = $j;
        }

        $mat{$i}{0} = $i;
    }

    # Some char-by-char processing is ahead, so prepare
    # array of chars from the strings
    #
    my @ar1 = split(//, $s1);
    my @ar2 = split(//, $s2);

    for (my $i = 1; $i <= $len1; ++$i)
    {
        for (my $j = 1; $j <= $len2; ++$j)
        {
            # Set the cost to 1 iff the ith char of $s1
            # equals the jth of $s2
            # 
            # Denotes a substitution cost. When the char are equal
            # there is no need to substitute, so the cost is 0
            #
            my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;

            # Cell $mat{$i}{$j} equals the minimum of:
            #
            # - The cell immediately above plus 1
            # - The cell immediately to the left plus 1
            # - The cell diagonally above and to the left + the cost
            #
            # We can either insert a new char, delete a char of
            # substitute an existing char (with an associated cost)
            #
            $mat{$i}{$j} = min_ref([$mat{$i-1}{$j} + 1,
                                $mat{$i}{$j-1} + 1,
                                $mat{$i-1}{$j-1} + $cost]);
        }
    }

    # Finally, the distance equals the rightmost bottom cell
    # of the matrix
    #
    # Note that $mat{$x}{$y} denotes the distance between the 
    # substrings 1..$x and 1..$y
    #
    return $mat{$len1}{$len2};
}


# minimal element of a list
#
sub min_ref
{
    my @list = @{$_[0]};
    my $min = $list[0];

    foreach my $i (@list)
    {
        $min = $i if ($i < $min);
    }

    return $min;
}


sub check_undef ($;$) {
  my $val = shift;
  my $replace = shift;
  $replace = defined $replace ? $replace : "";
  $val = defined $val ? $val : $replace;
  return $val;
}


sub reverseComplement ($) {
	my $seq = shift;
	(my $rc = reverse($seq)) =~ tr/ACGTacgtNn/TGCAtgcaNn/;
	return $rc;
}

sub listFilesInDir ($)
{
	my $path = shift;
	croak "Error: path $path not found" unless -d $path;
	my ($record, @data, @filenames);
	my @list = `ls -l $path`;

	foreach $record (@list) {
		next if $record =~ /^d/;
		@data = split(/\s+/,$record);
		next unless defined($data[8]);
		push(@filenames, join("/",$path,$data[8]));
	}
	
	carp "Warning: nothing found in path $path" unless @filenames > 0;

	return @filenames;
}

sub parseFilename ($) {
	my $fullname = shift;
	my ($name, $path, $ext) = fileparse($fullname, qr/\.\w{2,5}$/);
	return ($path, $name, $ext);
}

sub readChromsizeFile ($) {

	my $chrsizefile = shift;
	my %chrsize;

	# Open Chrsize file to get the total length of each Chr
	open (CHRSZ, "<", $chrsizefile) or croak "Error: cannot open $chrsizefile for reading";
	while (<CHRSZ>) {
		chomp;
		my ($chrnum, $len) = split(/\t/);
		next unless $chrnum =~ /[Cc]hr/;
		$chrsize{$chrnum} = $len;
	}

	return %chrsize;
}

#invert exit status of system call
sub System ($;$) {
  my $cmd = shift;
  my $quiet = shift;
  $quiet = 0 unless defined $quiet;
  print "$cmd\n" unless $quiet;
  my $status = system($cmd);
  return !$status;
}

sub Capture ($;$) {
  my $cmd = shift;
  my $quiet = shift;
  $quiet = 0 unless defined $quiet;
  print "$cmd\n" unless $quiet;
  my @output = capture($cmd);
  return @output;
}


sub dateFileFormat
{
  my ($sec, $min, $hour, $day,$month,$year) = localtime time;
  return( sprintf( "%04d%02d%02d_%02d%02d%02d",$year+1900,$month+1,$day,$hour,$min,$sec) );
}

sub cgiUpload ($$$) {
	my $cgi = shift;
	my $param = shift;
	my $outdir = shift;
	my $safe_chars =  "A-Za-z0-9_.-";
	my $file = $cgi->param($param);
	my ($path,$name,$ext) = parseFilename($file);
	$file = $name.$ext;
	$file =~ s/ /_/g;
	$file =~ s/[^$safe_chars]//g;
	$file = $outdir.$file;
	print $cgi->p("Uploading goodfile to server at $file...");

	my $fh  = $cgi->upload($param);
	if (defined $fh) {
		# Upgrade the handle to one compatible with IO::Handle:
		my $io_handle = $fh->handle;
		open OUTFILE,">",$file or croak "Error: could not open $file for writing";
		while (my $bytesread = $io_handle->read(my $buffer,1024)) {
			print OUTFILE $buffer;
		}
	}
	System("perl -pi -e 's/\r/\n/g' $file");

	return $file;
}

sub log2 ($) {
	my $n = shift;
	return log($n)/log(2);
}

sub log10 ($) {
	my $n = shift;
	return log($n)/log(10);
}

1;
