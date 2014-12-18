#!/usr/bin/perl

use Getopt::Long;
use POSIX qw(floor);

my $threshold = 10000;
my $info = 0;
my $numinf = 0;
my $means = 0;

my $resolution = 0;

my $code_dir = @ENV{'CODEDIR'};
my $data_dir = @ENV{'DATADIR'};

my $alpha_scale = 1.;
my $beta_scale = 1.;
my $scale = 1.;

my $no_header = 0;

GetOptions("threshold=i" => \$threshold,
	   "no-header" => \$no_header,
	   "info" => \$info,
	   "vinf" => \$num_inf,
	   "alpha-scale=f" => \$alpha_scale,
	   "beta-scale=f" => \$beta_scale,
	   "scale=f" => \$scale,
	   "means" => \$means,
	   "resolution=f" => \$resolution);


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

if (!($scale == 1.)) {
  $alpha_scale = $scale;
  $beta_scale = $scale;
}

my %data;
my %columns;

open(IN, $filename) or die "Can't read $filename\n";
while ($line = <IN>) {
  if ($line =~ /^[\w\+\-]+=([0-9\.]+) [\w\+\-]+=([0-9\.]+)/) {
    my $var1 = $1 / $alpha_scale;
    my $var2 = $2 / $beta_scale;
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
      if ($count > 0) {
        $data{$var2}{$var1} = $count;
      } else {
        $data{$var2}{$var1} = "-1";
      }
    }
    if ($resolution == 0 || abs($var1 / $resolution - floor($var1 / $resolution + 0.5)) < 1e-10) {
      $columns{$var1} = 1;
    }
  }
}

if ($info) {
  foreach $var2 (sort {$a <=> $b} keys %data) {
    foreach $var1 (sort keys %{$data{$var2}}) {
      print $var1."$delimiter".$data{$var2}{$var1}."\n";
    }
  }
} else {
  my $firstline = !$no_header;
  foreach $var2 (sort {$a <=> $b} keys %data) {
    if ($firstline) {
      $linestring = "   0"."$delimiter";
      foreach $var1 (sort {$a <=> $b} keys %columns) {
	if ($var1 < 10) {
	  $linestring .= sprintf("%.2f",$var1)."$delimiter";
	} elsif ($var1 < 100) {
	  $linestring .= sprintf("%.1f",$var1)."$delimiter";
	} else {
	  $linestring .= sprintf("%4.0f",$var1)."$delimiter";
	}
      }
      chop $linestring;
      print "$linestring\n";
      $firstline = 0;
    }

    $linestring = sprintf("%.2f",$var2)."$delimiter";
    foreach $var1 (sort {$a <=> $b} keys %columns) {
      $linestring .= sprintf("%4.0f",$data{$var2}{$var1})."$delimiter";
    }
    chop $linestring;
    print "$linestring\n";
  }
}

