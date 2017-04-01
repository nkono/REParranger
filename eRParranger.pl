####################################
# eRParranger ver 0.1
# 2017.3.31

use strict;
use Getopt::Long;

my %opts = ( limlen => 50000, output => "eRPoutput");
GetOptions(\%opts, qw( contig=s sam=s limlen=i output=s help) ) or exit 1;

if ($opts{help}) {
    die &show_help();
}

my $dieflag;
foreach my $field ( qw( contig sam ) ){
    $dieflag .= "Error: [$field] file is required.\n" if ! exists $opts{$field};
}
$dieflag .= "Error: [$opts{output}] directory already exists.\n" if -e $opts{output};
die $dieflag if $dieflag;

system("mkdir $opts{output}");

print "contig\t".$opts{contig}."\n";
print "sam\t".$opts{sam}."\n";
print "lim\t".$opts{limlen}."\n";

#####################################
######### Contigs from SPAdes #######
my %contigs;
my $key;
open(CONTIG, "$opts{contig}") or die "Could not open [$opts{contig}]\n";
while(<CONTIG>){
    chomp;
    if(/>(.+)/){
	$key = $1;
    }else{
	$contigs{$key}{SEQ} .= $_ if $key;
    }
}
close CONTIG;
die "Could not treat [$opts{contig}] as the FASTA format file\n" if scalar(keys %contigs) < 1;

#####################################
######### Coverage from SAM #########
my %sam;
my $unit = 1000;
open(SAM, $opts{sam}) or die "Could not open [$opts{sam}]\n";
while(<SAM>){
    chomp;
    my @data = split(/\t/, $_);
    my $contig = $data[2];
    my $len = length($contigs{$contig}{SEQ});
    next if $len < $opts{limlen};
    my $pos = int($data[3]/$unit) * $unit;
    $sam{$contig}{L} = $len;
    $sam{$contig}{$pos} ++;
}
close SAM;
die "Could not treat [$opts{sam}] as the SAM format file\n" if scalar(keys %sam) < 1;

my %coverage;
my %appres;
my $contig0;
my ($tpos, $maxcov) = (0, 0);

my $outappres = $opts{output}."/appres.tsv";
open(APRES, "> $outappres");
for my $contig (sort {$sam{$b}{L} <=> $sam{$a}{L}} keys %sam){
    for my $pos (sort {$a <=> $b} keys %{$sam{$contig}}){
	next if $pos eq "L";
	push(@{$coverage{$contig}{pos}}, $pos);
	push(@{$coverage{$contig}{cov}}, $sam{$contig}{$pos});
    }
    my ($slope, $intercept) = &approximation(\@{$coverage{$contig}{pos}}, \@{$coverage{$contig}{cov}});
    my ($pos1, $val1) = ($tpos, (0 * $slope + $intercept));
    $tpos += $sam{$contig}{L};
    my ($pos2, $val2) = ($tpos, ($sam{$contig}{L} * $slope + $intercept));
    $tpos ++;
    print APRES $contig."\t".$pos1."\t".$val1."\n";
    print APRES $contig."\t".$pos2."\t".$val2."\n";

    if($val1 < $val2){ # up
	$contigs{$contig}{STYLE} = "up";
	$appres{$contig}{Spos} = $pos2;
	$appres{$contig}{Epos} = $pos1;
	$appres{$contig}{Sval} = $val2;
	$appres{$contig}{Eval} = $val1;
    }elsif($val2 < $val1){ # down
	$contigs{$contig}{STYLE} = "down";
	$appres{$contig}{Spos} = $pos1;
	$appres{$contig}{Epos} = $pos2;
	$appres{$contig}{Sval} = $val1;
	$appres{$contig}{Eval} = $val2;
    }
#    print $contig."\t".$appres{$contig}{Sval}."\n";
    if($maxcov < $appres{$contig}{Sval}){
	$maxcov = $appres{$contig}{Sval};
	$contig0 = $contig;
    }
}
close APRES;

#####################################
#########     eRParrange     ########
my $flag = 0;
my %used;
my @next;

my @order;
$order[0] = $contig0;
while(1){
    my ($Sval0, $Eval0) = ($appres{$contig0}{Sval}, $appres{$contig0}{Eval});
    $used{$contig0} ++;
    my @contigs;
    if(scalar(@next) > 1){
	for(@next){
	    next if $used{$_};
	    push(@contigs, $_);
	}
    }else{
	for(keys %appres){
	    next if $used{$_};
	    push(@contigs, $_);
	}
    }
    my %regions;
    for my $contig (@contigs){
	my ($Sval1, $Eval1) = ($appres{$contig}{Sval}, $appres{$contig}{Eval});
	my $region;
	if($Eval0 < $Eval1){
	    next;
	}elsif($Sval0 > $Sval1 && $Sval1 > $Eval0){ # cross
	    $region = $Sval1 - $Eval0 + 1;
	}elsif($Eval0 > $Sval1){ # under
	    $region = 0;
	}
	$regions{$contig} = $region;
    }    
    @next = ();
    for my $contig (sort {$regions{$a} <=> $regions{$b}} keys %regions){
	if($regions{$contig} == 0){
	    push(@next, $contig);
	}else{
	    push(@next, $contig);
	    last;
	}
    }

    my $min_distance = 0;
    my $contig2;
    for my $contig (@next){
	if(!$min_distance){
	    $min_distance = abs($Eval0 - $appres{$contig}{Sval});
	    $contig2 = $contig;
	}elsif($min_distance > abs($Eval0 - $appres{$contig}{Sval})){
	    $min_distance = abs($Eval0 - $appres{$contig}{Sval});
	    $contig2 = $contig;
	}
    }
    last if !$contig2 && $flag;    
    if(!$contig2){
	$contig0 .= ":SWITCH:";
	$flag = 1;
	push(@order, "SWITCH");
	for my $contig (sort {$appres{$a}{Eval} <=> $appres{$b}{Eval}} keys %appres){
	    next if $used{$contig};
	    ($appres{$contig}{Sval}, $appres{$contig}{Eval}) = ($appres{$contig}{Eval}*(-1), $appres{$contig}{Sval}*(-1));
	}
    }else{
	push(@order, $contig2);
	$used{$contig2} ++;
	$contig0 .= ":".$contig2;
	if(scalar(@next) == 1){
	    if(!$flag){
		$contig0 .= ":SWITCH:";
		$flag = 1;
		push(@order, "SWITCH");
		for my $contig (sort {$appres{$a}{Eval} <=> $appres{$b}{Eval}} keys %appres){
		    next if $used{$contig};
		    ($appres{$contig}{Sval}, $appres{$contig}{Eval}) = ($appres{$contig}{Eval}*(-1), $appres{$contig}{Sval}*(-1));
		}
	    }
	}
    }
    $appres{$contig0}{Sval} = $Sval0;
    $appres{$contig0}{Eval} = $appres{$contig2}{Eval};
    $used{$contig0} ++;
}


my $outfasta = $opts{output}."/eRParranger.fasta";
my $outlist = $opts{output}."/eRParranger.list";
open(FASTA, "> $outfasta");
open(ORD, "> $outlist");
print FASTA ">tmp\n";
my $flag;
for my $contig (@order){
    if($contig eq "SWITCH"){
	$flag = 1;
	next;
    }
    if(!$flag){ # 0 # down
	if($contigs{$contig}{STYLE} eq "up"){
	    print FASTA complement($contigs{$contig}{SEQ})."\n";
	    print ORD $contig."\tup:down = comp\n";
	}elsif($contigs{$contig}{STYLE} eq "down"){
	    print FASTA $contigs{$contig}{SEQ}."\n";
	    print ORD $contig."\tdown:down = direct\n";
	}
    }else{      # 1 # up	
	if($contigs{$contig}{STYLE} eq "down"){
	    print FASTA complement($contigs{$contig}{SEQ})."\n";
	    print ORD $contig."\tdown:up = comp\n";
	}elsif($contigs{$contig}{STYLE} eq "up"){
	    print FASTA $contigs{$contig}{SEQ}."\n";
	    print ORD $contig."\tup:up = direct\n";
	}
    }
}
close FASTA;
close ORD;

sub approximation{
    my ($x, $y) = @_;
    my @x = @$x;
    my @y = @$y;
    my $num = scalar(@x);
    my ($sum_xy, $sum_xx, $sum_y, $sum_x);
    for my $i (0..$num-1){
        $sum_xy += $x[$i] * $y[$i];
        $sum_xx += $x[$i] * $x[$i];
        $sum_y += $y[$i];
        $sum_x += $x[$i];
    }
    my $ave_y = $sum_y / $num;
    my $ave_x = $sum_x / $num;
    my $slope = ($sum_xy - $num * $ave_y * $ave_x) / ($sum_xx - $num * $ave_x * $ave_x);
    my $intercept = $ave_y - $slope * $ave_x;
    return($slope, $intercept);
}


sub complement {
    my $nuc = reverse(shift);
    $nuc =~ tr
        [acgturymkdhbvwsnACGTURYMKDHBVWSN]
        [tgcaayrkmhdvbwsnTGCAAYRKMHDVBWSN];
    return $nuc;
}

sub show_help {
    my $help_doc = <<EOF;

Program: eRParranger (re-order system based on eRP curve)
Version: 0.1

Usage:   perl eRParranger <command> [options]

Command: -c FILE       SPAdes contig file (format: FASTA)
         -s FILE       read mapped file on SPAdes contig (format: SAM)
         -l INT        minimum contig length [50000]
	 -o STR        output dir name [eRPoutput]

Note: This program does not require other modules. On the other hand, users should prepare 
      two files (assembled contig file and read mapping file). The assembled contig file can be 
      obtained from SPAdes program (http://bioinf.spbau.ru/spades). The sequence reads mapped 
      SAM file on SPAdes contig can be generated by BWA program (http://bio-bwa.sourceforge.net).

EOF
return $help_doc;
}
