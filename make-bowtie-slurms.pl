#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <fastq-dir> <alignment-dir> <working-dir>\n" unless @ARGV==3;
my ($fastqDir,$alignmentDir,$workingDir)=@ARGV;

system("mkdir bowtie-slurms") unless -e "bowtie-slurms";
system("rm -f bowtie-slurms/*");

# Sa51-293T-15_S3_L001__FWD_paired.fq.gz
open(CMD,">bowtie-commands.sh");
my @files=`ls $fastqDir/*FWD_paired.fq.gz`;
my $index=1;
foreach my $file1 (@files) {
  chomp $file1;
  #next unless $file1=~/.*R1.*fastq/;
  $file1=~/([^\/]+)__FWD_paired.fq.gz/ || die $file1;
  my $ID=$1;
  my $file2="$fastqDir/${ID}__REV_paired.fq.gz";
  my $cmd="bowtie2 --local --maxins 12000 -N 0 --mp 100 --rdg 100,100 --rfg 100,100 -x bowtie/exon51 -1 $file1 -2 $file2 > $alignmentDir/$ID.sam";

  my $slurm="bowtie-slurms/$ID.slurm";
  open(OUT,">$slurm") || die $slurm;
  print OUT "#!/bin/bash
#
#SBATCH -J BOWTIE$index
#SBATCH -o BOWTIE$index.output
#SBATCH -e BOWTIE$index.output
#SBATCH -A BOWTIE$index
#
cd $workingDir
module load bowtie2/2.2.4-fasrc01
$cmd
";
  close(OUT);
  print CMD "$cmd\n";
  ++$index;
}
close(CMD);

open(OUT,">bowtie-slurms/submit.sh") || die;
print OUT "#!/bin/tcsh
ls *.slurm | perl -ne 'print \"sbatch \$_\"' | /bin/tcsh
";
close(OUT);


