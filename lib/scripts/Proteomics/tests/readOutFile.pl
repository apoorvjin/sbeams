#!/usr/local/bin/perl

###############################################################################
# Program     : readOutFile.pl
# Author      : Eric Deutsch <edeutsch@systemsbiology.org>
# $Id$
#
# Description : Test program to read a sequest .out file using
#               SBEAMS::Proteomics::Utilities
#
###############################################################################


  use strict;

  use lib qw (../perl ../../perl);
  use SBEAMS::Proteomics::Utilities;
  my $sbeamsPR = SBEAMS::Proteomics::Utilities->new();


  my ($i,$element,$key,$value);


  my $verbose = 0;
  my $inputfile = shift;

  unless ($inputfile) {
    if (1 == 1) {
      $inputfile =
        "/net/dblocal/data/macrogenics/data/CTCL/CTCL1/human_nci/".
        "CTCL1_0910_R01_042202/CTCL1_0910_R01_042202.3654.3654.2.out";
    } else {
      $inputfile =
        "/net/db/projects/proteomics/data/priska/ICAT/".
        "raftapr/raftapr_human/raft0052/raft0052.1052.1052.3.out";
    }
  }


  my %search_data = $sbeamsPR->readOutFile(inputfile => "$inputfile",
    verbose => "$verbose");


  print "\n\nsearch_data:\n";

  while ( ($key,$value) = each %search_data ) {
    printf("%22s = %s\n",$key,$value);
  }


  print "\nparameters:\n";
  while ( ($key,$value) = each %{$search_data{parameters}} ) {
    printf("%22s = %s\n",$key,$value);
  }


  print "\nmatches:\n";
  foreach $element ( @{$search_data{matches}} ) {
    print "  $element -> {reference}: $element->{reference}\n";
    if ($element->{search_hit_proteins}) {
      print "\t\t",join(',',@{$element->{search_hit_proteins}}),"\n";
    }
  }


