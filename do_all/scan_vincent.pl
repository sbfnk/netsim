#!/usr/bin/perl

use Getopt::Long;
use File::Path;
use File::Copy;

my $low_param = .0;
my $high_param = .0;
my $low_beta = .0;
my $high_beta = .0;
my $step = .02;

my $param_step;
my $beta_step;

my $param = "alpha";

my $code_dir = @ENV{'CODEDIR'};
my $data_dir = @ENV{'DATADIR'};

my $command = "$code_dir/do_all/do_stats.sh";

GetOptions("step=f" => \$step,
	   "param=s" => \$param,
	   "alpha-step=f" => \$param_step,
	   "beta-step=f" => \$beta_step);

my $usage = "Usage: scan_vincent.pl low_beta high_beta ".
  "low_param high_param [-s step] [-p param] [-a param-step] [-b beta-step]\n";

(@ARGV >= 4) || die ($usage);

$lb = shift(@ARGV);
$hb = shift(@ARGV);
$la = shift(@ARGV);
$ha = shift(@ARGV);

($lb <= $hb) && (la <= hb) || die ($usage); 

($param_step > 0) || ($param_step = $step);
($beta_step > 0) || ($beta_step = $step);
($param_step > 0) && ($beta_step > 0) || die($usage);

print "Varying $param in [$la,$ha) in steps of $param_step and beta in [$lb,$hb) in steps of $beta_step\n";
print "Parameters: @ARGV\n\n";

for ($alpha = $la; $alpha <= $ha; $alpha += $param_step) {
  for ($beta = $lb; $beta <= $hb; $beta += $beta_step) {
    print("$param=$alpha beta=$beta\n");
#    $arguments = sprintf("--$param=%.2f --beta=%.2f",$alpha,$beta);
    $arguments = "--$param=$alpha --beta=$beta";
    system("$command $arguments @ARGV");
  }
}

