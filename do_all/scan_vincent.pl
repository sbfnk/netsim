#!/usr/bin/perl

use Getopt::Long;
use File::Path;
use File::Copy;

my $low_param = .0;
my $high_param = .0;
my $low_beta = .0;
my $high_beta = .0;
my $step = .02;
my $cutoff = 0;

my $param_step;
my $beta_step;

my $param = "alpha";

my $code_dir = @ENV{'CODEDIR'};
my $data_dir = @ENV{'DATADIR'};

my $command = "$code_dir/do_all/do_stats.sh";

GetOptions("step=f" => \$step,
	   "param=s" => \$param,
	   "alpha-step=f" => \$param_step,
	   "beta-step=f" => \$beta_step,
	   "cutoff=f" => \$cutoff);

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

print "Varying $param in [$la,$ha) in steps of $param_step and beta in [$lb,$hb) in steps of $beta_step";
($cutoff > 0) && print ", cutoff $cutoff\n";
print "\n";
print "Parameters: @ARGV\n\n";

my $tmpfile = "/tmp/outbreaks.$$";

for ($alpha = $la; $alpha <= $ha; $alpha += $param_step) {
  my $davg = 0;
  for ($beta = $lb; $beta <= $hb && ($cutoff == 0 || $current_avg < $cutoff); $beta += $beta_step) {
#    print("$param=$alpha beta--=$beta beta-+=$beta\n");
#    $arguments = sprintf("--$param=%.2f --beta=%.2f",$alpha,$beta);
    $arguments = "--$param=$alpha --beta--=$beta --beta-+=$beta --r0 2";
    system("$command $arguments @ARGV > $tmpfile 2>&1");
    my $r0sum=0;
    my $r0count;
    my $dsum=0;
    my $dsqsum=0;
    my $dcount=0;
    my $isum=0;
    my $isqsum=0;
    my $icount=0;
    open (IN, "<$tmpfile");
    while (<IN>) {
      if (/^R0 generation 1 ([0-9\.]+) ([0-9]+)$/) {
        $r0sum += $2*$1;
        $r0count += $2;
      } elsif (/^Cumulative number of infections: ([0-9]+)$/) {
        $dsum += $1;
        $dsqsum += $1*$1;
        $dcount++;
      } elsif (/^Cumulative number of informations: ([0-9]+)$/) {
        $isum += $1;
        $isqsum += $1*$1;
        $icount++;
      }
    }
    close(IN);
    system("/bin/rm $tmpfile");
    if ($r0count == 0) { $r0count = 1; }
    if ($dcount == 0) { $dcount = 1; }
    if ($icount == 0) { $icount = 1; }
    $paramstr = sprintf("%.5f %.5f",$alpha,$beta);
    $davg = $dsum / ($dcount + 0.);
    $dvar = sqrt($dsqsum/($dcount + 0.) - $davg*$davg);
    $dstring = sprintf("%.0f (%.1f)",$davg,$dvar);
    $iavg = $isum / ($icount + 0.);
    $ivar = sqrt($isqsum/($icount + 0.) - $iavg*$iavg);
    $istring = sprintf("%.0f (%.1f)",$iavg,$ivar);
    $r0avg = $r0sum / ($r0count + 0.);
    $current_avg = $r0avg;
    $r0string = sprintf("%.3f",$r0avg);
    print "$paramstr $dstring $istring $r0string\n";
  }
}
