#!/usr/bin/perl

use Getopt::Long;
use File::Path;
use File::Copy;

my $lower = .0;
my $upper = .0;
my $step = .0;
my $param = "";
my $file_id = "";
my $force = 0;

my $code_dir = @ENV{'CODEDIR'};
my $data_dir = @ENV{'DATADIR'};

GetOptions("lower=f" => \$lower,
	   "upper=f" => \$upper,
	   "step=f" => \$step,
	   "params=s" => \@params,
	   "id=s" => \$file_id,
	   "force" => \$force);

my $usage = "Usage: vary_parameter.pl -f file_id -p param -l lower_bound ".
  "-u upper_bound -s step\n";

(($file_id ne "") && (@params > 0)) || die ($usage);
($code_dir ne "") || die("CODEDIR not set\n");
($data_dir ne "") || die("DATADIR not set\n");

@params = split(/,/,join(',',@params));

$file_id = "$file_id";

if ( -d "$data_dir/$file_id" ) {
  if ($force == 0) {
    die("$data_dir/$file_id exists. Use -f to overwrite\n");
  } else {
    rmtree("$data_dir/$file_id");
  }
}

mkdir("$data_dir/$file_id") or die ("Could not create $data_dir/$file_dir\n");

my %outfiles = ("mf" => "$data_dir/$file_id/$file_id.mf.dat",
#		"pa" => "$data_dir/$file_id/$file_id.pa.dat",
		"sim" => "$data_dir/$file_id/$file_id.sim.dat");

for ($i = $lower; $i < $upper; $i = $i + $step) {
  my $param_str = "";
  foreach $param (@params) {
    $param_str = "$param_str --$param=$i";
  }
  system("$code_dir/do_all/do_all.sh $file_id"."_$i -q $param_str".
	 " @ARGV");
  my %infiles = ("mf"=>"$data_dir/$file_id"."_$i/$file_id"."_$i.mf.dat",
#		 "pa"=>"$data_dir/$file_id"."_$i/$file_id"."_$i.pa.dat",
		 "sim"=>"$data_dir/$file_id"."_$i/$file_id"."_$i.sim.dat");
  while (($type, $filename) = each(%infiles)) {
    open(IN, $filename) or die "Can't read $filename\n";
    my $saveline = "";
    while ($line = <IN>) {
      if ($line =~ /\S/) {
	$saveline = $line;
      }
    }
    close(IN);

    @data = split(/\s/, $saveline);
    $data[0] = $i;

    open(OUT, "+>>$outfiles{$type}") or die "Can't write to $outfiles{$type}";
    print OUT "@data\n";
    close(OUT);

  }

}

copy("$data_dir/$file_id"."_$lower/"."$file_id"."_$lower.gp",
     "$data_dir/$file_id/$file_id.gp") or die("Could not copy gnuplot file");

system("$code_dir/do_all/makevary.bash $file_id $param");
