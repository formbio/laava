#!/usr/bin/perl -w
#merge_megadepthcts.pl

my $prefix = shift @ARGV;
my $ctfile = shift @ARGV;

my %gene;
open BED, "<$ctfile" or die $!;
while (my $line = <BED>) {
    chomp($line);
    my ($chr,$start,$end,$name,$numreads,$numbases,$length,$fracov) = split(/\t/,$line);
    my $key = $chr.':'.$start.'-'.$end;
    my @names = split(/,/, $name);
    my %uni;
    foreach $exon (@names) {
	my ($gene,$ensembl,$trxid,$exonnum,$strand) = split(/\|/,$exon);
	next unless ($ensembl);
	$uni{$ensembl} = $gene;
    }
    foreach $ens (keys %uni) {
	$sym = $uni{$ens};
	$genect{join("|",$ens,$sym)} += $numreads;
    }
}

open OUT, ">$prefix\.bedtools.cov.txt" or die $!;
foreach $gene (keys %genect) {
    my ($ens,$sym) = split(/\|/,$gene);
    print OUT join("\t",$ens,$sym,$genect{$gene}),"\n";
}
