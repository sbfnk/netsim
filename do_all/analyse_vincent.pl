#!/usr/bin/perl

use Getopt::Long;

my $threshold = 10000;
my $info = 0;
my $numinf = 0;
my $means = 0;

my $code_dir = @ENV{'CODEDIR'};
my $data_dir = @ENV{'DATADIR'};

my $no_header = 0;

GetOptions("threshold=i" => \$threshold,
	   "no-header" => \$no_header,
	   "info" => \$info,
	   "vinf" => \$num_inf,
	   "means" => \$means);


my $usage = "Usage: analyse_vincent.pl [-t threshold] [-x x-axis] [-i] [-v]".
            " [-m] filename\n";

$filename = shift(@ARGV);

($filename ne "") || die($usage);

my $linestring = "";

my $delimiter = " ";
my $countstring;

if ($info) {
  $countstring = "informations";
} elsif ($num_inf) {
  $countstring = "infected vertices";
} else {
  $countstring = "infections";
}

my %data;

open(IN, $filename) or die "Can't read $filename\n";
while ($line = <IN>) {
  if ($line =~ /^[\w\+\-]+=([0-9\.]+) [\w\+\-]+=([0-9\.]+)/) {
    my $var1 = $1;
    my $var2 = $2;
    my $count = 0;
    my $sum = 0;
    $line = <IN>;
    while ($line =~ /\S/) {
      if ($line =~ /$countstring.*\D(\d+)\w*$/) {
	if ($means) {
	  $sum += $1;
	  $count++;
	} elsif ($1 >= $threshold) { $count++; }
      }
      $line = <IN>;
    }
    if ($means) {
      $data{$var2}{$var1} = ($sum / $count);
    } else {
      $data{$var2}{$var1} = $count;
    }
  }
}

if ($info) {
  foreach $var2 (sort keys %data) {
    foreach $var1 (sort keys %{$data{$var2}}) {
      print $var1."$delimiter".$data{$var2}{$var1}."\n";
    }
  }
} else {
  my $firstline = !$no_header;
  foreach $var2 (sort keys %data) {
    if ($firstline) {
      $linestring = "   0"."$delimiter";
      foreach $var1 (sort keys %{$data{$var2}}) {
	$linestring .= sprintf("%.2f",$var1)."$delimiter";
      }
      chop $linestring;
      print "$linestring\n";
      $firstline = 0;
    }

    $linestring = sprintf("%.2f",$var2)."$delimiter";
    foreach $var1 (sort keys %{$data{$var2}}) {
      $linestring .= sprintf("%4.0f",$data{$var2}{$var1})."$delimiter";
    }
    chop $linestring;
    print "$linestring\n";
  }
}

