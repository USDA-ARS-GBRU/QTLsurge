#!/usr/bin/env perl
#
# Author: Brian Abernathy
# Date: 2019-12-09

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $vcf_file;
my $output_prefix;
my $pop1_name;
my $pop2_name;
my $par1_name;
my $par2_name;
my $pop1_field;
my $pop2_field;
my $par1_field;
my $par2_field;
my $avg_pop1_depth;
my $avg_pop2_depth;
my $avg_par1_depth;
my $avg_par2_depth;
my $min_pop1_depth;
my $min_pop2_depth;
my $min_par1_depth;
my $min_par2_depth;
my $max_pop1_depth;
my $max_pop2_depth;
my $max_par1_depth;
my $max_par2_depth;
my $depth_factor = 2;
my $min_depth = 10;
my $max_depth;
my $min_qual = 50;
my $min_mq = 50;
my $pop_ratio_filter;
my $help;

&parse_args();
&get_vcf_fields();

$output_prefix =~ s/\.$//;

my $output_file = "$output_prefix.vcf";
my $stats_file = "$output_prefix.stats";

open(OUTPUT, '>', $output_file) || &gen_error("cannot write $output_file: $!");
open(STATS, '>', $stats_file) || &gen_error("cannot write $stats_file: $!");

my $low_qual = 0;
my $low_mq = 0;
my $low_depth = 0;
my $high_depth = 0;
my $contam = 0;
my $parents_het = 0;
my $parents_not_poly = 0;
my $high_pop_ratio = 0;
my $passed = 0;

open(VCF, '<', $vcf_file) or &gen_error("cannot read $vcf_file: $!");

while (my $line = <VCF>) {
	chomp($line);

	if ($line =~ /^#/) {
		if ($line =~ /^#CHROM/) {
			print(OUTPUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$pop1_name\t$pop2_name");

			if (defined($par1_name)) {
				print(OUTPUT "\t$par1_name");
			}

			if (defined($par2_name)) {
				print(OUTPUT "\t$par2_name");
			}

			print(OUTPUT "\n");
		}

		else {
			print(OUTPUT "$line\n");
		}

		next();
	}

	my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $gt_format, @vcf_fields) = split(/\t/, $line);

	if ($qual < $min_qual) {
		$low_qual++;

		next();
	}


	my $mq;

	if ($info =~ /MQ\=(\d+)\;/) {
		$mq = $1;
	}

	if (! defined($mq) || $mq < $min_mq) {
		$low_mq++;

		next();
	}


	my $pop1_vcf = $vcf_fields[$pop1_field];
	my $pop2_vcf = $vcf_fields[$pop2_field];

	my ($pop1_depth, $pop1_ad, $pop1_gt) = get_depth_alleles($pop1_vcf, $gt_format);
	my ($pop2_depth, $pop2_ad, $pop2_gt) = get_depth_alleles($pop2_vcf, $gt_format);

	my $par1_vcf;
	my $par2_vcf;
	
	my ($par1_depth, $par1_ad, $par1_gt);
	my ($par2_depth, $par2_ad, $par2_gt);

	if (defined($par1_field)) {
		$par1_vcf = $vcf_fields[$par1_field];

		($par1_depth, $par1_ad, $par1_gt) = get_depth_alleles($par1_vcf, $gt_format);
	}
	
	if (defined($par2_field)) {
		$par2_vcf = $vcf_fields[$par2_field];

		($par2_depth, $par2_ad, $par2_gt) = get_depth_alleles($par2_vcf, $gt_format);
	}


	if ($pop1_depth < $min_depth || $pop2_depth < $min_depth || (defined($par1_depth) && $par1_depth < $min_depth) || (defined($par2_depth) && $par2_depth < $min_depth)) {
		$low_depth++;

		next();
	}

	if ((defined($min_pop1_depth) && $pop1_depth < $min_pop1_depth) || (defined($min_pop2_depth) && $pop2_depth < $min_pop2_depth)) {
		$low_depth++;

		next();
	}

	if (defined($par1_depth) && defined($min_par1_depth) && $par1_depth < $min_par1_depth) {
		$low_depth++;

		next();
	}

	if (defined($par2_depth) && defined($min_par2_depth) && $par2_depth < $min_par2_depth) {
		$low_depth++;

		next();
	}

	if (defined($max_depth)) {
		if ($pop1_depth > $max_depth || $pop2_depth > $max_depth) {
			$high_depth++;

			next();
		}

		if ((defined($par1_depth) && $par1_depth > $max_depth) || (defined($par2_depth) && $par2_depth > $max_depth)) {
			$high_depth++;

			next();
		}
	}

	if ((defined($max_pop1_depth) && $pop1_depth > $max_pop1_depth) || (defined($max_pop2_depth) && $pop2_depth > $max_pop2_depth)) {
		$high_depth++;

		next();
	}

	if (defined($par1_depth) && defined($max_par1_depth) && $par1_depth > $max_par1_depth) {
		$high_depth++;

		next();
	}

	if (defined($par2_depth) && defined($max_par2_depth) && $par2_depth > $max_par2_depth) {
		$high_depth++;

		next();
	}


	if (defined($par1_gt) && defined($par2_gt) && ($par1_gt eq '0/0' && $par2_gt eq '0/0') && ($pop1_gt ne '0/0' || $pop2_gt ne '0/0')) {
		$contam++;
	}

	if (defined($par1_gt) && defined($par2_gt)) {
		my $par1_gt_homo = 0;
		my $par2_gt_homo = 0;

		my ($par1_gt1, $par1_gt2) = split(/\//, $par1_gt);

		if ($par1_gt1 eq $par1_gt2) {
			$par1_gt_homo = 1;
		}

		my ($par2_gt1, $par2_gt2) = split(/\//, $par2_gt);

		if ($par2_gt1 eq $par2_gt2) {
			$par2_gt_homo = 1;
		}


    	if ($par1_gt_homo == 0 || $par2_gt_homo == 0) {
			$parents_het++;

			next();
		}
	

		if ($par1_gt eq $par2_gt) {
			$parents_not_poly++;

			next();
		}
	}


	my $ratio_fail = 0;

	if (defined($pop_ratio_filter)) {
		my @pop1_ads = split(',', $pop1_ad);
		my @pop2_ads = split(',', $pop2_ad);

		foreach my $index (0..$#pop1_ads) {
			if ($ratio_fail == 1) {
				next();
			}

			my $pop1_ratio = $pop1_ads[$index] / $pop1_depth * 100;
			my $pop2_ratio = $pop2_ads[$index] / $pop2_depth * 100;

			if ($pop1_ratio > 95 && $pop2_ratio > 95) {
				$high_pop_ratio++;
				$ratio_fail = 1;

				next();
			}
		}
	}

	if ($ratio_fail == 1) {
		next();
	}

	
	print(OUTPUT join("\t", $chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $gt_format, $pop1_vcf, $pop2_vcf));

	if (defined($par1_vcf)) {
		print(OUTPUT "\t$par1_vcf");
	}

	if (defined($par2_vcf)) {
		print(OUTPUT "\t$par2_vcf");
	}

	print(OUTPUT "\n");
	
	$passed++;
}


my $total = $low_qual + $low_mq + $low_depth + $high_depth + $parents_het + $parents_not_poly + $high_pop_ratio + $passed;

if ($total > 0) {
	print(STATS "low qual: $low_qual\t", sprintf("%.2f", $low_qual / $total * 100), "%\n");
	print(STATS "low mq: $low_mq\t", sprintf("%.2f", $low_mq / $total * 100), "%\n");
	print(STATS "low depth: $low_depth\t", sprintf("%.2f", $low_depth / $total * 100), "%\n");
	print(STATS "high depth: $high_depth\t", sprintf("%.2f", $high_depth / $total * 100), "%\n");
	print(STATS "both parents het: $parents_het\t", sprintf("%.2f", $parents_het / $total * 100), "%\n");
	print(STATS "parents not poly: $parents_not_poly\t", sprintf("%.2f", $parents_not_poly / $total * 100), "%\n");
	print(STATS "high pop ratio: $high_pop_ratio\t", sprintf("%.2f", $high_pop_ratio / $total * 100), "%\n");
	print(STATS "possible contam: $contam\t", sprintf("%.2f", $contam / $total * 100), "%\n");
	print(STATS "passed: $passed \t", sprintf("%.2f", $passed / $total * 100), "%\n");
	print(STATS "total: $total\n");
}

close(OUTPUT);
close(STATS);

exit(0);


sub get_depth_alleles {
	my $sample = shift();
	my $format = shift();

	my $sample_depth;

	my @format_fields = split(':', $format);
	my @sample_fields = split(':', $sample);

	my $gt;
	my $ad;

	foreach my $field_num (0..$#format_fields) {
		my $format_val = $format_fields[$field_num];

		if ($format_val eq 'GT') {
			$gt = $sample_fields[$field_num];
		}

		if ($format_val eq 'AD') {
			$ad = $sample_fields[$field_num];
		}
	}

	if (! defined($gt)) {
		&gen_error("VCF record does not contain GT field");
	}

	if (! defined($ad)) {
		&gen_error("VCF record does not contain AD field");
	}

	my @allele_depths = split(',', $ad);

	foreach my $allele_index (0..$#allele_depths) {
		my $allele_depth = $allele_depths[$allele_index];

		$sample_depth += $allele_depth;
	}

	return($sample_depth, $ad, $gt);
}


sub get_vcf_fields {
	my %vcf_fields = ();

	open(VCF, '<', $vcf_file) or &gen_error("cannot read $vcf_file: $!");

	while (my $line = <VCF>) {
		if ($line =~ /^#CHROM/) {
			chomp($line);

			my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $gt_format, @vcf_fields) = split(/\t/, $line);

			foreach my $field_num (0..$#vcf_fields) {
				my $field_val = $vcf_fields[$field_num];

				$vcf_fields{$field_val} = $field_num;
			}

			last();
		}
	}

	close(VCF);

	if (defined($pop1_name)) {
		if (exists($vcf_fields{$pop1_name})) {
			$pop1_field = $vcf_fields{$pop1_name};
		}

		else {
			&gen_error("population 1 name: $pop1_name not found in VCF header");
		}
	}

	if (defined($pop2_name)) {
		if (exists($vcf_fields{$pop2_name})) {
			$pop2_field = $vcf_fields{$pop2_name};
		}

		else {
			&gen_error("population 2 name: $pop2_name not found in VCF header");
		}
	}

	if (defined($par1_name)) {
		if (exists($vcf_fields{$par1_name})) {
			$par1_field = $vcf_fields{$par1_name};
		}

		else {
			&gen_error("parent 1 name: $par1_name not found in VCF header");
		}
	}

	if (defined($par2_name)) {
		if (exists($vcf_fields{$par2_name})) {
			$par2_field = $vcf_fields{$par2_name};
		}

		else {
			&gen_error("parent 2 name: $par2_name not found in VCF header");
		}
	}

	return(0);
}


sub parse_args {
	my $result = GetOptions ('v|vcf=s' => \$vcf_file,
							'o|output=s' => \$output_prefix,
							'pop1_name=s' => \$pop1_name,
							'pop2_name=s' => \$pop2_name,
							'par1_name=s' => \$par1_name,
							'par2_name=s' => \$par2_name,
							'pop1_depth=f' => \$avg_pop1_depth,
							'pop2_depth=f' => \$avg_pop2_depth,
							'par1_depth=f' => \$avg_par1_depth,
							'par2_depth=f' => \$avg_par2_depth,
							'min_depth=f' => \$min_depth,
							'max_depth=f' => \$max_depth,
							'q|qual=f' => \$min_qual,
							'mq=f' => \$min_mq,
							'depth_factor=f' => \$depth_factor,
							'pop_ratio' => \$pop_ratio_filter,
							'h|help' => \$help);

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 99);
	}

	if (! defined($vcf_file)) {
		&arg_error('vcf file not specified');
	}

	if (! defined($output_prefix)) {
		&arg_error('output prefix not specified');
	}

	if (! defined($pop1_name)) {
		&arg_error('population 1 name not specified');
	}

	if (! defined($pop2_name)) {
		&arg_error('population 2 name not specified');
	}

	if (defined($avg_pop1_depth) || defined($avg_pop2_depth) || defined($avg_par1_depth) || defined($avg_par2_depth)) {
		if (! defined($depth_factor)) {
			&arg_error('depth factor not specified');
		}

		if ($depth_factor <= 0) {
			&arg_error('depth factor must be > 0');
		}

		if (defined($avg_pop1_depth)) {
			if ($avg_pop1_depth <= 0) {
				&arg_error('population 1 depth must be > 0');
			}

			$min_pop1_depth = $avg_pop1_depth / $depth_factor;
			$max_pop1_depth = $avg_pop1_depth * $depth_factor;
		}

		if (defined($avg_pop2_depth)) {
			if ($avg_pop2_depth <= 0) {
				&arg_error('population 2 depth must be > 0');
			}

			$min_pop2_depth = $avg_pop2_depth / $depth_factor;
			$max_pop2_depth = $avg_pop2_depth * $depth_factor;
		}

		if (defined($avg_par1_depth)) {
			if ($avg_par1_depth <= 0) {
				&arg_error('parent 1 depth must be > 0');
			}

			$min_par1_depth = $avg_par1_depth / $depth_factor;
			$max_par1_depth = $avg_par1_depth * $depth_factor;
		}

		if (defined($avg_par2_depth)) {
			if ($avg_par2_depth <= 0) {
				&arg_error('parent 2 depth must be > 0');
			}

			$min_par2_depth = $avg_par2_depth / $depth_factor;
			$max_par2_depth = $avg_par2_depth * $depth_factor;
		}
	}

	if (defined($min_depth) && $min_depth < 1) {
		&arg_error('minimum depth must be >= 1');
	}

	if (defined($max_depth) && $max_depth <= 0) {
		&arg_error('maximum depth must be > 0');
	}

	if (! defined($min_qual)) {
		&arg_error('minimum qual not specified');
	}

	if ($min_qual < 0) {
		&arg_error('minimum qual must be >= 0');
	}

	if (! defined($min_mq)) {
		&arg_error('minimum mq not specified');
	}

	if ($min_mq < 0) {
		&arg_error('minimum mq must be >= 0');
	}

	return(0);
}


sub arg_error {
	my $msg = shift();

	if (defined($msg)) {
		print("error: $msg\n\n");
	}

	pod2usage(-exitbal => 0, -verbose => 99);
}


sub gen_error {
	my $msg = shift();

	if (defined($msg)) {
		print("error: $msg\n");
	}

	exit(0);
}


__END__

=head1 NAME

qtl.filter.vcf.pl

=head1 SYNOPSIS

qtl.filter.vcf.pl [options]

=head1 DESCRIPTION

qtl.filter.vcf.pl filters a VCF file for further QTL processing. The VCF file must contain data for 2 populations and optionally 2 parents.

=head1 OPTIONS

 -v --vcf        (required) input vcf file

 -o --output     (required) output file prefix, used as the base filename for ouput vcf and stats files

 --pop1_name     (required) population 1 name as defined in VCF header

 --pop2_name     (required) population 2 name as defined in VCF header

 --par1_name     parent 1 name as defined in VCF header

 --par2_name     parent 2 name as defined in VCF header

 --min_depth     minimum depth, supercedes other depth filters (default: 10)

 --max_depth     maximum depth, supercedes other depth filters (default: no max)

 -q --qual       minimum SNP QUAL value (default: 50)

 --mq            minimum SNP MQ value (default: 50)

 --pop_ratio     enable filter removing SNPs where both bulks are >95% for the same allele (possibly as a result of segregation distortion)

 -h --help       display help menu

The options below set dynamic min/max depth filters for each average population/parent depth specified. These filters are not enabled by default.

min depth = average depth / depth factor

max depth = average depth * depth factor

 --pop1_depth    population 1 average depth

 --pop2_depth    population 2 average depth

 --par1_depth    parent 1 average depth

 --par2_depth    parent 2 average depth

 --depth_factor  value to divide/multiply average depths by to determine min/max depths (default: 2);

=cut
