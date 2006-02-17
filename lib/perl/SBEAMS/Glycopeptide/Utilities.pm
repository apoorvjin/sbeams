package SBEAMS::Glycopeptide::Utilities;

sub new {
  my $class = shift;
  my $this = {};
  bless $this, $class;
  return $this;
}

sub match_count {
  my $self = shift;
  my %args = @_;
  print STDERR "pep: $args{pepseq}\n";
  return unless $args{pepseq} && $args{protseq};

  my @cnt = split( $args{pepseq}, $args{protseq}, -1 );
  return $#cnt;
}

# map_peptide_to_protein, 0-based
sub map_peptide_to_protein {
	my $self = shift;
	my %args = @_;
	my $pep_seq = $args{pepseq};
	my $protein_seq = $args{protseq};
	
	if ( $protein_seq =~ /$pep_seq/ ) {
		my $start_pos = length($`);    
		my $stop_pos = length($pep_seq) + $start_pos;  
		return ($start_pos, $stop_pos);	
	}else{
		return;
	}
}

sub getDigestPeptide {
  my $self = shift;
  my %args = @_;
  for my $req ( qw( begin end protseq ) ) {
    die "Missing required parameter $req" unless defined $args{$req};
  }
  my $length =  $args{end} - $args{begin};
  my $seq = '';
  if ( !$args{begin} ) {
    $seq = substr( '-' . $args{protseq}, $args{begin}, $length + 2 );
  } elsif ( $args{end} == length($args{protseq}) ) {
    $seq = substr( $args{protseq} . '-' , $args{begin} -1, $length + 2 );
  } else {
    $seq = substr( $args{protseq}, $args{begin} -1, $length + 2 );
  }
  $seq =~ s/^(.)(.*)(.)$/$1\.$2\.$3/;
  return $seq;
}

# Returns reference to an array holding the 0-based indices of a pattern 'X'
# in the peptide sequence
sub get_site_positions {
  my $self = shift;
  my %args = @_;
  $args{pattern} ||= 'N.[S|T]';
  return unless $args{seq};

  my @posn;
  while ( $args{seq} =~ m/$args{pattern}/g ) {
    my $posn = length($`);
    push @posn, $posn;# pos($string); # - length($&) # start position of match
  }
  print STDERR "Found $posn[0] for NxS/T in $args{seq}\n";
  return \@posn;
}

sub get_current_prophet_cutoff {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my $cutoff = $sbeams->getSessionAttribute( key => 'glyco_prophet_cutoff' );
  $cutoff = 0.5 if !defined $cutoff;
  return $cutoff;
}

sub set_prophet_cutoff {
  my $self = shift;
  my $cutoff = shift || return;
  my $sbeams = $self->getSBEAMS();
  $sbeams->setSessionAttribute( key => 'glyco_prophet_cutoff',
                              value => $cutoff );
  return 1;
}

sub clean_pepseq {
  my $this = shift;
  my $seq = shift || return;
  $seq =~ s/\-MET\(O\)/m/g;
  $seq =~ s/N\*/n/g;
  $seq =~ s/N\#/n/g;
  $seq =~ s/M\#/m/g;
  $seq =~ s/d/n/g;
  $seq =~ s/U/n/g;

  # Trim of leading/lagging amino acids
  $seq =~ s/^.\.//g;
  $seq =~ s/\..$//g;
  return $seq;
}

###############################################################################
# getResidueMasses: Get a hash of masses for each of the residues
###############################################################################
sub getResidueMasses {
  my %args = @_;
  my $SUB_NAME = 'getResidueMasses';

  #### Define the residue masses
  my %residue_masses = (
    I => 113.1594,   # Isoleucine
    V =>  99.1326,   # Valine
    L => 113.1594,   # Leucine
    F => 147.1766,   # Phenyalanine
    C => 103.1388,   # Cysteine
    M => 131.1926,   # Methionine
    A =>  71.0788,   # Alanine
    G =>  57.0519,   # Glycine
    T => 101.1051,   # Threonine
    W => 186.2132,   # Tryptophan
    S =>  87.0782,   # Serine
    Y => 163.1760,   # Tyrosine
    P =>  97.1167,   # Proline
    H => 137.1411,   # Histidine
    E => 129.1155,   # Glutamic_Acid (Glutamate)
    Q => 128.1307,   # Glutamine
    D => 115.0886,   # Aspartic_Acid (Aspartate)
    N => 114.1038,   # Asparagine
    K => 128.1741,   # Lysine
    R => 156.1875,   # Arginine

    X => 113.1594,   # L or I
    B => 114.5962,   # avg N and D
    Z => 128.6231,   # avg Q and E
    U => 100,        # ?????

  );

  return \%residue_masses;

}

sub getPeptideMass {
  my $self = shift;
  my %args = @_;
  die "Missing required parameter sequence" unless $args{sequence};
  $self->{_rmass} ||= $self->getResidueMasses();
  my $seq = uc( $args{sequence} );
  my @seq = split( "", $seq );
  my $mass = 0;
  foreach my $r ( @seq ) {
    if ( !defined $self->{_rmass}->{$r} ) {
      print STDERR "Undefined residue $r is getPeptideMass\n";
      $self->{_rmass}->{$r} = $self->{_rmass}->{U} # Assign 'average' mass.
    }
    $mass += $self->{_rmass}->{$r};
  }
  print STDERR "$args{sequence} => $mass\n";
  return $mass;
}


1;

