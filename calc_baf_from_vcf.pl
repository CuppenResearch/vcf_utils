#! /usr/bin/perl

use strict;

my $snpsift = "/hpc/cog_bioinf/common_scripts/snpEff_v4/SnpSift.jar";
my $vcf = shift;

my %filehandles;

my $header = `grep "^#CHROM" -m 1 $vcf`;
chomp($header);
my @columns = split/\s+|\t/, $header;
my @samples = @columns[9..$#columns];
foreach my $sample (@samples) {
  open $filehandles{$sample}, ">", "$sample\_baf_vcf.txt";
}

open F, "java -jar -Xmx4G $snpsift extractFields $vcf CHROM POS REF ALT GEN[*].AD | ";
while(<F>) {
  chomp;
  next if $_ =~ /^#/;
  my ($chr, $pos, $ref, $alt, @ad) = split/\t/, $_;
  foreach my $index (0..$#ad) {
    my ($ref_depth, $alt_depth) = ($1, $2) if $ad[$index] =~ /(\d+),(\d+)/;
    my $baf = "NaN";
    $baf = ($alt_depth/($ref_depth+$alt_depth)) if ($ref_depth+$alt_depth) > 0;
    my $sample = $samples[$index];
    print {$filehandles{$sample}} join("\t", $chr, $pos, $baf) . "\n";
  }
}
close F;

foreach my $sample (keys %filehandles) {
  close $filehandles{$sample};
  system "Rscript -e 'data <- read.table(\"$sample\_baf_vcf.txt\", sep=\"\\t\"); pdf(\"$sample\_baf_vcf.pdf\",width=20,height=10); for (chr in c(1:22,\"X\",\"Y\")) { if (nrow(data[data[,1] == chr,]) > 0) { plot(data[data[,1] == chr,2:3],ylim=c(0,1),pch=20,main=chr,ylab=\"BAF\",xlab=\"Position\"); }; }; dev.off()'\n";
}


