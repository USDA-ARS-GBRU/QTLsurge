#!/usr/bin/perl -w
use strict;
use warnings;
my $VERSION="0.1";
use Data::Dumper;

#use Math::Random;

my $db = shift; #vcf file
my $cycle = shift;
open IN, "$db";

# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	yellow.2018.w22.sort.bam	purple.2018.w22.sort.bam
# CM004348.2	114	.	C	G	28.3555	.	.	GT:PL:AD	0/0:0,9,131:3,0	0/1:61,0,126:6,4
# CM004348.2	10720	.	A	G	999	.	.	GT:PL:AD	1/1:255,39,0:0,13	1/1:255,69,0:0,23


print "chr\tposition\thighBulkFreq\tlowBulkFreq\tdeltaSNP\tcycle\n";
my $ids;
my $oldChrom = "";
while (my $l = <IN>) {
	next if $l =~ /^##/;
	chomp $l;
	if ($l =~ /^#CHROM/) {
		my @s = split /\t/, $l;
		my @fields = splice(@s,0,9);
		unshift(@s, "snp");
		$ids = join("\t", @s)."\n";
		next;
	}
	my @s = split /\t/, $l;
	my $chrom = $s[0];
	my $pos = $s[1];
	my $ref = $s[3];
	my @alts = split /\,/, $s[4];
	my $format = $s[8];
	die("Your genotype format is not \"GT:PL:AD\".  Looks like you used a different pipeline to generate the vcf file? You can modifiy this script if possible to support your format.  Otherwise please follow the pipline given at https://github.com/USDA-ARS-GBRU/QTLsurge") unless $format eq "GT:PL:AD";
	my @fields = splice(@s,0,9);
	my @toPrint = ($chrom,$pos);
	my $zero = 0;
	foreach (@s) {
		my @s2 = split /\:/, $_;
		my @s3 = split /\,/, $s2[2];
		#skips lines without sum
		if ($s3[0] + $s3[1] == 0) {
			$zero = 1;
			last;
		}
		my $freq = $s3[0] / ($s3[0] + $s3[1]);
		push(@toPrint, $freq);
	}
	next if $zero;
	push(@toPrint, $toPrint[3] - $toPrint[2]);
	push(@toPrint, $cycle);
	print(join("\t", @toPrint)."\n");
}
close IN;



