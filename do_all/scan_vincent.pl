#!/usr/bin/perl

use Getopt::Long;
use File::Path;
use File::Copy;

my $low_alpha = .0;
my $high_alpha = .0;
my $low_beta = .0;
my $high_beta = .0;
my $step = .02;

my $code_dir = @ENV{'CODEDIR'};
my $data_dir = @ENV{'DATADIR'};

GetOptions("step=f" => \$step);

my $usage = "Usage: vary_parameter.pl low_beta high_beta ".
  "low_alpha high_alpha [-s step]\n";

(@ARGV >= 4) || die ($usage);

$lb = shift(@ARGV);
$hb = shift(@ARGV);
$la = shift(@ARGV);
$ha = shift(@ARGV);

($lb <= $hb) && (la <= hb) || die ($usage);

for ($alpha = $la; $alpha <= $ha; $alpha += $step) {
  for ($beta = $lb; $beta <= $hb; $beta += $step) {
    print "alpha=$alpha beta=$beta\n";
    $arguments = "--alpha=$alpha --beta-+=$beta --beta--=$beta ".
      "@ARGV";
    $command = "$code_dir/do_all/do_stats.sh";
    system "$command $arguments";
  }
}

