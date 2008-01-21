#!/usr/bin/perl

use Getopt::Long;

my $no_header = 0;

my $usage = "Usage: info_dist.pl filename\n";

$filename = shift(@ARGV);

($filename ne "") || die($usage);

my $delimiter = " ";

my $sum = .0;
my $sqsum = .0;
my $count = 0;

open(IN, $filename) or die "Can't read $filename\n";
while ($line = <IN>) {
  if ($line =~ /[0-9\.\-e]+[^0-9\.\-e]+([0-9\.\-e]+)/) {
    my $number = $1;
    $sum += $number;
    $sqsum += $number*$number;
    $count++;
  }
}

my $mean;
my $dev;

if ($count > 0) {
  $mean = $sum / $count;
  $dev = sqrt(($sqsum - 2*$sum*$mean + $count*$mean*$mean) / $count);
} else {
  $mean = 0;
  $dev = 0;
}

print "$mean $dev\n";
