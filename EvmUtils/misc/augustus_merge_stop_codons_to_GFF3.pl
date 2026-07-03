#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $model_type = "Augustus";

my $usage = "usage: $0 augustus.gff.output\n\n";

my $input_file = $ARGV[0] or die $usage;

main: {
	# per model: list of [end5, end3, feat_type] segments, feat_type is 'CDS' or 'stop_codon'
	my %segments;

	## parse input file
	open (my $fh, $input_file) or die "Error, cannot open file $input_file";
	while (<$fh>) {
		if (/^\#/) { next; }
		chomp;
		unless (/\w/) { next; }

		my @x = split(/\t/);
		if ($x[2] eq 'CDS' || $x[2] eq "stop_codon") {
			my $scaffold = $x[0];
			my $orient = $x[6];
			my $lend = $x[3];
			my $rend = $x[4];

			my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

			my $info = $x[8];
			$info =~ /Parent=([^;\s]+)/;
			my $model_id = $1 or die "Error, cannot parse model from $info";

            #print "MODEL: [$model_id]\n";

			my $model = join("$;", $scaffold, $model_id);

			push (@{$segments{$model}}, [$end5, $end3, $x[2]]);
		}
	}

	close $fh;


	## Generate gff3 output
	foreach my $model (keys %segments) {

		my ($scaffold, $model_id) = split(/$;/, $model);

		my $coords_href = merge_stop_codon_into_CDS($segments{$model});

		my $gene_obj = new Gene_obj();
		$gene_obj->populate_gene_object($coords_href, $coords_href);
		$gene_obj->{asmbl_id} = $scaffold;
		$gene_obj->{TU_feat_name} = "gene.$model_id";
		$gene_obj->{Model_feat_name} = "model.$model_id";
		$gene_obj->{com_name} = "$model_type prediction";

        $gene_obj->join_adjacent_exons();

		print $gene_obj->to_GFF3_format(source => $model_type) . "\n";

        # print $gene_obj->toString() . "\n\n";

	}


	exit(0);

}


####
# Augustus can be run with --stopCodonExcludedFromCDS=True (stop_codon is a
# separate segment immediately adjacent to the CDS; join_adjacent_exons()
# later merges it in) or --stopCodonExcludedFromCDS=False (the stop codon is
# already included within the CDS segment, and the stop_codon feature is
# purely redundant annotation of those same 3 bp). A stop_codon segment must
# only be kept when it doesn't already overlap a CDS segment for this model;
# otherwise it introduces a spurious overlapping/backwards CDS exon.
sub merge_stop_codon_into_CDS {
    my ($segs_aref) = @_;

    my @cds_segs   = grep { $_->[2] eq 'CDS' } @$segs_aref;
    my @stop_segs  = grep { $_->[2] eq 'stop_codon' } @$segs_aref;

    my %coords;
    foreach my $seg (@cds_segs) {
        my ($end5, $end3) = @$seg;
        $coords{$end5} = $end3;
    }

    foreach my $seg (@stop_segs) {
        my ($end5, $end3) = @$seg;
        my ($slo, $shi) = sort { $a <=> $b } ($end5, $end3);

        my $overlaps_existing_cds = 0;
        foreach my $cds_seg (@cds_segs) {
            my ($clo, $chi) = sort { $a <=> $b } @$cds_seg[0,1];
            if ($slo <= $chi && $clo <= $shi) {
                $overlaps_existing_cds = 1;
                last;
            }
        }

        next if $overlaps_existing_cds;

        $coords{$end5} = $end3;
    }

    return \%coords;
}
		
		
