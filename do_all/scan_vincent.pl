#!/usr/bin/perl

use Getopt::Long;
use File::Path;
use File::Copy;

my $low_alpha = .0;
my $high_alpha = .0;
my $low_beta = .0;
my $high_beta = .0;
my $step = .02;

my $alpha_step;
my $beta_step;

my $code_dir = @ENV{'CODEDIR'};
my $data_dir = @ENV{'DATADIR'};

my $command = "$code_dir/do_all/do_stats.sh";

GetOptions("step=f" => \$step,
	   "alpha-step=f" => \$alpha_step,
	   "beta-step=f" => \$beta);

my $usage = "Usage: vary_parameter.pl low_beta high_beta ".
  "low_alpha high_alpha [-s step] [-a alpha-step] [-b beta-step]\n";

(@ARGV >= 4) || die ($usage);

$lb = shift(@ARGV);
$hb = shift(@ARGV);
$la = shift(@ARGV);
$ha = shift(@ARGV);

($lb <= $hb) && (la <= hb) || die ($usage); 

($alpha_step > 0) || ($alpha_step = $step);
($beta_step > 0) || ($beta_step = $step);
($alpha_step > 0) && ($beta_step > 0) || die($usage);

print "Varying alpha in [$la,$ha] and beta in [$lb,$hb]\n";
print "Parameters: @ARGV\n\n";

for ($alpha = $la; $alpha <= $ha; $alpha += $step) {
  for ($beta = $lb; $beta <= $hb; $beta += $step) {
    printf("alpha=%.2f beta=%.2f\n", $alpha, $beta);
    $arguments = sprintf("--alpha=%.2f --beta-+=%.2f --beta--=%.2f",$alpha,$beta,$beta);
    system("$command $arguments @ARGV");
  }
}

