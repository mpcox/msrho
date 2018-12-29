#!/usr/bin/env perl

# Programme: msrho v.2.0
# A programme to calculate Forster's rho from ms simulations
# Institution: Arizona Research Laboratories, University of Arizona
# Author: Murray Cox
# Date: May 2007 (v.1.0: 21 November 2006)
# Developed on perl v5.8.6 built for darwin-thread-multi-2level

# Checked that code still runs
# Date: December 2018
# Confirmed on perl v5.18.2 built for darwin-thread-multi-2level

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil floor);
use FileHandle;

# set up command line options
my $density_flag;
GetOptions(
	"d=s" => \$density_flag		# If defined, hold file name for rho(i) density output
);

# usage information
my $usage = "error: correct usage is: ms ... -T | msrho [-d density_output_file]\n";

# make density file
if( defined $density_flag ){
	if ( -e $density_flag ) {
		die $usage;
	} else {
		open(DENSITY,">$density_flag");
	}
}

# disable buffering on STDOUT and DENSITY filehandles
STDOUT->autoflush(1);
if( defined $density_flag ){
	DENSITY->autoflush(1);
}

my $sim_size = 0;
my $seg_sites;
my $chr_c = 0;

my @matrix;
my @identity;
my @mutations;

my $line_number = 0;

MAIN: while ( <STDIN> ) {
	
	# prepare line
	my $line = $_;
	chomp($line);
	
	# check for simulation size (first line)
	$sim_size = $1 if $line =~ /ms (\d+)/;
	die $usage if $sim_size == 0;
	
	# print header line
	if( $line_number == 0 ) {
		print STDOUT "S\trho\tvar\tsd\tl_95\tu_95\n";
	}

	# step to matrix
	if( $line =~ /\/\// ){
		undef($seg_sites);
		$chr_c = 0;
		next MAIN;
	}
	if( $line =~ /segsites: (\d+)/ ){
		$seg_sites = $1;

		# if there are no segregating sites
		if( $seg_sites == 0 ){
			print STDOUT "0\t0\t0\t0\t0\t0\n";
			undef($seg_sites);
		}
		
		next MAIN;
	}
	
	# extract matrix
	my @line_s = split(//, $line);
	my $complete_flag = 0;
	
	if(    defined($seg_sites)
			&& scalar(@line_s) == $seg_sites
	    && $seg_sites != 0
	    && ($line_s[0] == 0 || $line_s[0] == 1) ) {
		
		# check if identical to previous line in matrix
		for( my $i = 0; $i < scalar(@matrix); $i++) {
			if( $line =~ $matrix[$i] ){
				$identity[$i]++;
				$chr_c++;
				$complete_flag = 1;
				next MAIN if $chr_c != $sim_size;
			}
		}
		
		# else push onto matrix
		if( $complete_flag == 0 ) {
			push(@matrix, $line );
			push(@identity, 1 );
			$chr_c++;
			next MAIN if $chr_c != $sim_size;
		}
	}
	
	if( $chr_c == $sim_size ) {
		
		# calculate number of mutations
		foreach my $entry (@matrix) {
			push(@mutations, $entry =~ tr/1//);
		}
		
		# calculate rho statistic
		my $rho_mean_sum = 0;
		for( my $a = 0; $a < scalar(@identity); $a++) {
			$rho_mean_sum += $identity[$a] * $mutations[$a];
		}
		my $rho_mean = $rho_mean_sum / $sim_size;
		
		# calculate rho variance
		my $rho_var_sum = 0;
		for( my $a = 0; $a < scalar(@identity); $a++) {
			$rho_var_sum += ($identity[$a] ** 2) * $mutations[$a];
		}
		my $rho_var = $rho_var_sum / ($sim_size ** 2);
		
		# calculate rho standard deviation
		my $rho_sd = sqrt($rho_var);

		# calculate confidence intervals
		my @density;
		for( my $b = 0; $b < scalar(@identity); $b++) {
			for( my $c = 0; $c < $identity[$b]; $c++) {
				push(@density, $mutations[$b]);
				
				if( $density_flag ) {
					print DENSITY $mutations[$b], "\t";
				}
				
			}
		}
		
		if( $density_flag ) {
			print DENSITY "\n";
		}
		
		my @sorted_density = sort { $a <=> $b } @density;
		
		my $lower1 = floor(0.025 * scalar(@sorted_density));
		my $lower2 =  ceil(0.025 * scalar(@sorted_density));
		my $upper1 = floor(0.975 * scalar(@sorted_density)) - 1;
		my $upper2 =  ceil(0.975 * scalar(@sorted_density)) - 1;
		
	 	my $lower_CI  = ($sorted_density[$lower1] + $sorted_density[$lower2]) / 2;
  		my $upper_CI  = ($sorted_density[$upper1] + $sorted_density[$upper2]) / 2;
		
		# print results to screen
		print STDOUT "$seg_sites\t$rho_mean\t$rho_var\t$rho_sd\t$lower_CI\t$upper_CI\n";
		
		# reset variables
		$chr_c = 0;
		undef(@matrix);
		undef(@identity);
		undef(@mutations);
		undef(@density);
		
	}

	$line_number++;
}

# successful programme termination
exit 0;

