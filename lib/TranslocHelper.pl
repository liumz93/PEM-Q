# Source for larger subroutines having to do with major program flow

use Switch::Plain;

sub wrap_alignment ($$) {

  my $pe = shift;

  my $aln = shift;

  return undef unless defined $aln;
  croak "Error: first argument to wrap_alignment must be either \"R1\" or \"R2\"" unless ($pe eq "R1" || $pe eq "R2");
  # my %wrapper :shared;
  my $wrapper = {Qname => $aln->qname};

  if ($wrapper->{Qname} =~ /(\d+:\d+:\d+$)/ ) {
    $wrapper->{QnameShort} = $1;
  } else {
    $wrapper->{QnameShort} = $wrapper->{Qname};
  }

  $wrapper->{ID} = Data::GUID->new->as_string;

  if ($aln->unmapped) {
    $wrapper->{Unmapped} = 1;
    $wrapper->{Seq} = $pe eq "R1" ? $aln->query->dna : reverseComplement($aln->query->dna);
    $wrapper->{Qual} = $pe eq "R1" ? $aln->query->qscore : [reverse @{$aln->query->qscore}];
    $wrapper->{Qstart} = 0;
    return $wrapper;
  }

  $wrapper->{Rname} = $aln->seq_id;
  $wrapper->{Rstart} = $aln->start;
  $wrapper->{Rend} = $aln->end;
  $wrapper->{AS} = $aln->aux =~ /AS:i:(\d+)/ ? $1 : 0;
  $wrapper->{MDZ} = $aln->aux =~ /MD:Z:(\S+)/ ? $1 : 0;

  $wrapper->{Unmapped} = 0;

  if ($pe eq "R1") {
    $wrapper->{Strand} = $aln->strand;
  } else {
    $wrapper->{Strand} = -1 * $aln->strand;
  }


  if ($wrapper->{Strand} == 1) {
    $wrapper->{Seq} = $aln->query->dna;
    $wrapper->{Qual} = $aln->query->qscore;
    $wrapper->{Qstart} = $aln->query->start;
    $wrapper->{Qend} = $aln->query->end;
    # $wrapper->{CigarA} = $aln->cigar_array;
  } else {
    $wrapper->{Seq} = reverseComplement($aln->query->dna);
    $wrapper->{Qual} = [reverse @{$aln->query->qscore}];
    $wrapper->{Qstart} = $aln->l_qseq - $aln->query->end + 1;
    $wrapper->{Qend} = $aln->l_qseq - $aln->query->start + 1;
    # $wrapper->{CigarA} = [reverse @{$aln->cigar_array}];
  }

  $wrapper->{CigarA} = $wrapper->{MDZ} =~ /\d[ACGT]/ ?
                         remap_cigar($aln->cigar_array,$aln->query->dna,$aln->dna) :
                         $aln->cigar_array;

  $wrapper->{CigarA} = [reverse @{$wrapper->{CigarA}}] unless $wrapper->{Strand} == 1;

  $wrapper->{Cigar} = cigar_array_to_string($wrapper->{CigarA});
  $wrapper->{Qlen} = length($wrapper->{Seq});

  return $wrapper;

}



sub alignment_matches_simulation ($$) {
  my $aln = shift;
  my $sim_aln = shift;

  if (defined $aln->{Rname} &&
      $aln->{Rname} eq $sim_aln->{Rname} &&
      $aln->{Strand} == $sim_aln->{Strand}) {
    my $sim_start;
    my $sim_end;
    if ($aln->{Strand} == 1) {
      $sim_start = $sim_aln->{Junction};
      $sim_end = $sim_aln->{Junction} + $sim_aln->{PreyLen} - 1;
    } else {
      $sim_start = $sim_aln->{Junction} - $sim_aln->{PreyLen} + 1;
      $sim_end = $sim_aln->{Junction};
    }
    debug_print("comparing aln ".$aln->{Rstart}."-".$aln->{Rend}." to sim ".$sim_start."-".$sim_end,4,$aln->{QnameShort});
    return 1 if ($aln->{Rstart} <= $sim_end && $aln->{Rend} >= $sim_start);
  }
  return 0;
}



sub find_optimal_coverage_set ($$) {

  my $R1_alns_ref = shift;
  my $R2_alns_ref = shift;


  my @graph = ();
  my $OCS_ptr;

  my @R1_alns = sort {$a->{Qstart} <=> $b->{Qstart}} shuffle(values $R1_alns_ref);
  my @R2_alns = sort {$a->{Qstart} <=> $b->{Qstart}} shuffle(values $R2_alns_ref);

  my $qname = $R1_alns[0]->{QnameShort};

  debug_print("finding OCS",2,$qname);

  debug_print("searching through R1 and R1/R2 pairs",3,$qname);
  foreach my $R1_aln (@R1_alns) {

    next if $R1_aln->{Unmapped};

    # my $graphsize = scalar @graph;

    my $new_node = {R1 => $R1_aln};

    my $init_score = score_edge($new_node);
    $new_node->{score} = $init_score if defined $init_score;

    foreach my $node (@graph) {
      my $edge_score = score_edge($node,$new_node);
      next unless defined $edge_score;

      if (! exists $new_node->{score} || $edge_score > $new_node->{score}) {
        $new_node->{score} = $edge_score;
        $new_node->{back_ptr} = $node;

        debug_print("new best edge score for R1 node; set backpointer",3,$qname);
      }
    }

    if (defined $new_node->{score}) {
      push(@graph,$new_node);
      if (! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score}) {
        $OCS_ptr = $new_node;
        debug_print("new top OCS score for R1 node; set OCS pointer",3,$qname);
      }
    }

    foreach my $R2_aln (@R2_alns) {

      next unless defined $R2_aln && ! $R2_aln->{Unmapped};

      next unless pair_is_proper($R1_aln,$R2_aln);

      my $new_pe_node = {R1 => $R1_aln, R2 => $R2_aln};

      my $init_score = score_edge($new_pe_node);
      $new_pe_node->{score} = $init_score if defined $init_score;
      # if ($graphsize > 0) {
        # foreach my $node (@graph[0..($graphsize-1)]) {
        foreach my $node (@graph) {

          my $edge_score = score_edge($node,$new_pe_node);
          next unless defined $edge_score;

          if (! defined $new_pe_node->{score} || $edge_score > $new_pe_node->{score}) {
            $new_pe_node->{score} = $edge_score;
            $new_pe_node->{back_ptr} = $node;

            debug_print("new best edge score for R1/R2 node; set backpointer",3,$qname);

          }
        }
      # }
      if (defined $new_pe_node->{score}) {
        push(@graph,$new_pe_node);
        if (! defined $OCS_ptr || $new_pe_node->{score} > $OCS_ptr->{score}) {
          $OCS_ptr = $new_pe_node;
          debug_print("new top OCS score for R1/R2 node; set OCS pointer",3,$qname);
        }
      }
    }
  }


  debug_print("searching through R2 alignments",3,$qname);
  
  foreach my $R2_aln (@R2_alns) {

    next unless defined $R2_aln && ! $R2_aln->{Unmapped};

    my $new_node = {R2 => $R2_aln};


    foreach my $node (@graph) {
      my $edge_score = score_edge($node,$new_node);
      next unless defined $edge_score;

      if (! defined $new_node->{score} || $edge_score > $new_node->{score}) {
        $new_node->{score} = $edge_score;
        $new_node->{back_ptr} = $node;
        debug_print("new best edge score for R2 node; set backpointer",3,$qname);
      }
    }

    if (defined $new_node->{score}) {
      push(@graph,$new_node) ;
      if (! defined $OCS_ptr || $new_node->{score} > $OCS_ptr->{score}) {
        $OCS_ptr = $new_node;
        debug_print("new top OCS score for R2 node; set OCS pointer",3,$qname);
      }
    }

  }

  if (defined $OCS_ptr) {

    my @OCS = ($OCS_ptr);

    while (defined $OCS_ptr->{back_ptr}) {
      $OCS_ptr->{back_ptr}->{score} = $OCS_ptr->{score};
      $OCS_ptr = $OCS_ptr->{back_ptr};
      unshift(@OCS,$OCS_ptr)
    }
    debug_print("found OCS of length ".scalar @OCS,3,$qname);
    return \@OCS;

  } else {

    debug_print("no defined OCS pointer, returning unmapped OCS",3,$qname);
    my @unmapped_OCS = ( { R1 => { Qname => $R1_alns[0]->{Qname},
                                   QnameShort => $R1_alns[0]->{QnameShort},
                                   Seq => $R1_alns[0]->{Seq},
                                   Qual => $R1_alns[0]->{Qual},
                                   Unmapped => 1 },
                           R2 => { Qname => $R2_alns[0]->{Qname},
                                   QnameShort => $R2_alns[0]->{QnameShort},
                                   Seq => $R2_alns[0]->{Seq},
                                   Qual => $R2_alns[0]->{Qual},
                                   Unmapped => 1 } } );
    return \@unmapped_OCS;

  }

}

sub score_edge ($;$) {
  my $node1 = shift;
  my $node2 = shift;

  my $score;
  my $qname;

  if (defined $node1->{R1}) {
    $qname = $node1->{R1}->{QnameShort};
  } else {
    $qname = $node1->{R2}->{QnameShort};
  }

  debug_print("scoring edge",3,$qname);
  debug_print("node1:R1:".stringify_alignment($node1->{R1}),4,$qname)
    if defined $node1->{R1};
  debug_print("node1:R2:".stringify_alignment($node1->{R2}),4,$qname)
    if defined $node1->{R2};

  if (defined $node2) {

    debug_print("node2:R1:".stringify_alignment($node2->{R1}),4,$qname)
      if defined $node2->{R1};
    debug_print("node2:R2:".stringify_alignment($node2->{R2}),4,$qname)
      if defined $node2->{R2};


    my $Rname1;
    my $Rname2;
    my $query_intersection;
    my $query_union;

    if ( defined $node1->{R2} ) {

      if (defined $node2->{R1}) {
        debug_print("score=undef; node1 has R2, node2 must not have R1",4,$qname);
        return undef;
      }
      unless (defined $node2->{R2}) {
        debug_print("score=undef; node1 has R2, node2 must have R2",4,$qname);
        return undef;
      }
      if ($node1->{R2}->{Rname} eq "Adapter") {
        debug_print("score=undef; node1 is adapter, nothing may follow",4,$qname);
        return undef;
      }
      if ($node2->{R2}->{ID} eq $node1->{R2}->{ID}) {
        debug_print("score=undef; node1 and node2 have equivalent R2",4,$qname);
        return undef;
      }
      unless ($node2->{R2}->{Qstart} >= $node1->{R2}->{Qstart}) {
        debug_print("score=undef; node1 R2 Qstart greater than node2 R2 Qstart",4,$qname);
        return undef;
      }
      unless ($node2->{R2}->{Qend} >= $node1->{R2}->{Qend}) {
        debug_print("score=undef; node1 R2 Qend greater than node2 R2 Qend",4,$qname);
        return undef;
      }
      $Rname1 = $node1->{R2}->{Rname};
      $Rname2 = $node2->{R2}->{Rname};
      $query_intersection = aln_query_overlap($node1->{R2},$node2->{R2});
      $query_union = $node2->{R2}->{Qend} - $node1->{R2}->{Qstart} + 1;


    } else {

      unless (defined $node2->{R1}) {
        debug_print("score=undef; node1 has no R2, node2 must have R1",4,$qname);
        return undef;
      }
      if ($node1->{R1}->{Rname} eq "Adapter") {
        debug_print("score=undef; node1 is adapter, nothing may follow",4,$qname);
        return undef;
      }
      if ($node2->{R1}->{ID} eq $node1->{R1}->{ID}) {
        debug_print("score=undef; node1 and node2 have equivalent R1",4,$qname);
        return undef;
      }
      unless ($node2->{R1}->{Qstart} >= $node1->{R1}->{Qstart}) {
        debug_print("score=undef; node1 R1 Qstart greater than node2 R1 Qstart",4,$qname);
        return undef;
      }
      unless ($node2->{R1}->{Qend} >= $node1->{R1}->{Qend}) {
        debug_print("score=undef; node1 R1 Qend greater than node2 R1 Qend",4,$qname);
        return undef;
      }
      $Rname1 = $node1->{R1}->{Rname};
      $Rname2 = $node2->{R1}->{Rname};
      $query_intersection = aln_query_overlap($node1->{R1},$node2->{R1});
      $query_union = $node2->{R1}->{Qend} - $node1->{R1}->{Qstart} + 1;

    }

    # node1 and node2 pass tests, calc score

    my $overlap_correction;    
    my $brk_pen;

    if ($Rname2 eq "Adapter") {
      $overlap_correction = 0;
      $brk_pen = 0;
    } else {

      if ($query_intersection/$query_union > $main::params->{max_overlap}) {
        debug_print("score=undef; node 1 and node 2 overlap greater than maximum",4,$qname);
        return undef;
      }

      $overlap_correction = $main::params->{match_award} * $query_intersection;
      
      # Add in correction for brksite cassette?
      $brk_pen = $main::params->{brk_pen};
    }

    my $R1_AS = defined $node2->{R1} ? $node2->{R1}->{AS} : 0;
    my $R2_AS = defined $node2->{R2} ? $node2->{R2}->{AS} : 0;
    
    my $PE_pen = 0;

    if (defined $node2->{R1} && defined $node2->{R2}) {

      # Strange case:
      # If node2 R2 makes a proper pair with node1 R1,
      # But node2 R1 does not make a proper pair with node1 R1,
      # Do not accept
      if (pair_is_proper($node1->{R1},$node2->{R2})) {
        
        if (! pair_is_proper($node1->{R1},$node2->{R1})) {
          debug_print("score=undef; guessing node2 R2 is pair of node1 R1",4,$qname);
          return undef;
        }
      }


      $PE_pen = calculate_paired_end_penalty($node2->{R1},$node2->{R2});

    }

    debug_print("node1 score: ".$node1->{score},4,$qname);
    debug_print("node2 R1 AS: ".$R1_AS,4,$qname);
    debug_print("node2 R2 AS: ".$R2_AS,4,$qname);
    debug_print("node2 PE penalty: ".$PE_pen,4,$qname);
    debug_print("Junction penalty: ".$brk_pen,4,$qname);
    debug_print("Overlap correction: ".$overlap_correction,4,$qname);

    $score = $node1->{score} + $R1_AS + $R2_AS - $PE_pen - $brk_pen - $overlap_correction;
    
    debug_print("Edge score: ".$score,4,$qname);


  } else {
    # First node in OCS

    unless (defined $node1->{R1}) {
      debug_print("score=undef; first node must have R1",4,$qname);
      return undef;
    }

    if ($node1->{R1}->{Rname} eq "Adapter") {
      debug_print("score=undef; first node cannot be adapter",4,$qname);
      return undef;
    }

    if ($main::params->{force_bait}) {
      # Penalize for alignments away from bait site

      unless ($node1->{R1}->{Rname} eq $main::params->{brksite}->{aln_name}) {
        debug_print("score=undef; first node must be on brksite chr",4,$qname);
        return undef;
      }
      unless ($node1->{R1}->{Strand} == $main::params->{brksite}->{aln_strand}) {
        debug_print("score=undef; first node must be on brksite strand",4,$qname);
        return undef;
      }

      my $brk_start_dif = $main::params->{brksite}->{aln_strand} == 1 ?
                            abs($node1->{R1}->{Rstart} - $main::params->{brksite}->{primer_start}) :
                            abs($node1->{R1}->{Rend} - $main::params->{brksite}->{primer_start}) ;

      unless ($brk_start_dif < $main::params->{max_brkstart_dif}) {
        debug_print("score=undef; first node must be within max-brkstart-dif",4,$qname);
        return undef;
      }

    }

    # node1 passes tests; calc score

    my $R1_AS = $node1->{R1}->{AS};
    my $R2_AS = defined $node1->{R2} ? $node1->{R2}->{AS} : 0;

    my $PE_pen = 0;

    if (defined $node1->{R1} && defined $node1->{R2}) {
      $PE_pen = calculate_paired_end_penalty($node1->{R1},$node1->{R2});
    }

    debug_print("node1 R1 AS: ".$R1_AS,4,$qname);
    debug_print("node1 R2 AS: ".$R2_AS,4,$qname);
    debug_print("node1 PE penalty: ".$PE_pen,4,$qname);

    $score = $R1_AS + $R2_AS - $PE_pen;

    debug_print("Edge score: ".$score,4,$qname);

  }

  return $score;

}


sub calculate_paired_end_penalty ($$) {
  my $R1_aln = shift;
  my $R2_aln = shift;

  my $PE_gap = $R1_aln->{Strand} == 1 ?
                $R2_aln->{Rstart} - $R1_aln->{Rend} - 1 :
                $R1_aln->{Rstart} - $R2_aln->{Rend} - 1 ;
  
  my $PE_pen = 0;

  my $main_penalty = $PE_gap > 0 ? $main::params->{brk_pen} + $main::params->{pe_pen} * $PE_gap/$main::params->{max_pe_gap} : $main::params->{brk_pen};

  $PE_pen += $main_penalty;

  # Correct for overlap (half credit)
  # my $overlap_penalty = $main::params->{match_award} * aln_reference_overlap($R1_aln,$R2_aln) / 2;
  # $PE_pen += $overlap_penalty;

  # Penalize for unused R1
  my $unused_R1_penalty = $main::params->{match_award} * ($R1_aln->{Qlen} - $R1_aln->{Qend});
  $PE_pen += $unused_R1_penalty;

  debug_print("paired-end penalties: main:".$main_penalty." unused_R1:".$unused_R1_penalty,4,$R1_aln->{QnameShort});

  return $PE_pen;

}


sub create_tlxl_entries ($) {

  my $OCS_ref = shift;
  my @OCS = @$OCS_ref;

  debug_print("creating tlxl",2,$OCS[0]->{R1}->{QnameShort});

  my @tlxls = ();


  if ( $OCS[0]->{R1}->{Unmapped} ) {
    my $tlxl = { Qname => $OCS[0]->{R1}->{Qname},
                 QnameShort => $OCS[0]->{R1}->{QnameShort},
                 R1_Seq => $OCS[0]->{R1}->{Seq},
                 R2_Seq => $OCS[0]->{R2}->{Seq},
                 Unmapped => 1 };
    return [$tlxl];
  }


  foreach my $Qseg (@OCS[0..$#OCS]) {
    my $tlxl = {};

    if (defined $Qseg->{R1}) {
      $tlxl->{Qname} = $Qseg->{R1}->{Qname};
      $tlxl->{QnameShort} = $Qseg->{R1}->{QnameShort};
      $tlxl->{OCS_score} = $Qseg->{score};
      $tlxl->{R1_Qstart} = $Qseg->{R1}->{Qstart};
      $tlxl->{R1_Qend} = $Qseg->{R1}->{Qend};
      $tlxl->{R1_Qlen} = $Qseg->{R1}->{Qlen};
      $tlxl->{Rname} = $Qseg->{R1}->{Rname};
      $tlxl->{Strand} = $Qseg->{R1}->{Strand};
      $tlxl->{R1_Rstart} = $Qseg->{R1}->{Rstart};
      $tlxl->{R1_Rend} = $Qseg->{R1}->{Rend};
      $tlxl->{R1_Rgap} = $Qseg->{R1_Rgap};
      $tlxl->{R1_Cigar} = $Qseg->{R1}->{Cigar};
      $tlxl->{R1_CigarA} = $Qseg->{R1}->{CigarA};
      $tlxl->{R1_Seq} = $Qseg->{R1}->{Seq};
      $tlxl->{R1_Qual} = $Qseg->{R1}->{Qual};
      $tlxl->{R1_ID} = $Qseg->{R1}->{ID};
      $tlxl->{R1_AS} = $Qseg->{R1}->{AS};
    }

    if (defined $Qseg->{R2}) {

      unless (defined $Qseg->{R1}) {
        $tlxl->{Qname} = $Qseg->{R2}->{Qname};
        $tlxl->{QnameShort} = $Qseg->{R2}->{QnameShort};
        $tlxl->{OCS_score} = $Qseg->{score};
        $tlxl->{Rname} = $Qseg->{R2}->{Rname};
        $tlxl->{Strand} = $Qseg->{R2}->{Strand};        
      }
      $tlxl->{R2_Qstart} = $Qseg->{R2}->{Qstart};
      $tlxl->{R2_Qend} = $Qseg->{R2}->{Qend};
      $tlxl->{R2_Qlen} = $Qseg->{R2}->{Qlen};
      $tlxl->{R2_Rstart} = $Qseg->{R2}->{Rstart};
      $tlxl->{R2_Rend} = $Qseg->{R2}->{Rend};
      $tlxl->{R2_Rgap} = $Qseg->{R2_Rgap};
      $tlxl->{R2_Cigar} = $Qseg->{R2}->{Cigar};
      $tlxl->{R2_CigarA} = $Qseg->{R2}->{CigarA};
      $tlxl->{R2_Seq} = $Qseg->{R2}->{Seq};
      $tlxl->{R2_Qual} = $Qseg->{R2}->{Qual};
      $tlxl->{R2_ID} = $Qseg->{R2}->{ID};
      $tlxl->{R2_AS} = $Qseg->{R2}->{AS};

    }

    push(@tlxls,$tlxl);

  }

  return \@tlxls;

}


sub create_tlx_entries ($$) {

  my $tlxls = shift;
  my $refs = shift;

  debug_print("creating tlx",2,$tlxls->[0]->{QnameShort});

  my @tlxs;

  # Unmapped reads are easy
  # TLX will be one element long with simply the Qname and nothing else
  if ($tlxls->[0]->{Unmapped}) {
    my $tlx = {};
    $tlx->{Qname} = $tlxls->[0]->{Qname};
    $tlx->{QnameShort} = $tlxls->[0]->{QnameShort};
    $tlx->{JuncID} = 1;
    push(@tlxs,$tlx);
    return(\@tlxs);
  }

  # We need to merge the alignments if there is an R1/R2 match
  # And return a query coordinate map for R1 -> merged seq
  # and for R2 -> merged seq
  my $Qseq;
  my $Qlen;

  my $R1_map;
  my $R2_map;
  my $Cigar_array_ref;

  # First decide if we even need a merge
  my $R2_idx = firstidx {defined $_->{R2_Seq}} @$tlxls;

  if ($R2_idx < 0) {


    if ($tlxls->[$#$tlxls]->{Rname} eq "Adapter") {
      $Qseq = substr($tlxls->[0]->{R1_Seq},0,$tlxls->[$#$tlxls]->{R1_Qend});
      $Qlen = length($Qseq);
    } else {
      $Qseq = $tlxls->[0]->{R1_Seq};
      $Qlen = "";    
    }

    # Set R1 as Qseq and just a one-to-one map for R1 if no merge
    $R1_map = [1..length($Qseq)];
    $R2_map = [];
    $Cigar_array_ref = [];
  } else {
    # print "merging alignments on $R2_idx\n";
    ($Qseq,$R1_map,$R2_map,$Cigar_array_ref) = merge_alignments($tlxls,$refs) ;
    $Qlen = length($Qseq);
  }


  foreach my $i (0..$#$tlxls) {
    # print "TLX for ".$i."th segment\n";
    

    my $tlx = {};
    my $b_tlxl = $tlxls->[$i];

    $tlx->{Qname} = $b_tlxl->{Qname};
    $tlx->{QnameShort} = $b_tlxl->{QnameShort};
    $tlx->{B_Rname} = $b_tlxl->{Rname};
    $tlx->{B_Strand} = $b_tlxl->{Strand};
    $tlx->{Seq} = $Qseq;
    $tlx->{JuncID} = $i + 1;

    # 
    $tlx->{Qlen} = $Qlen;

  
    # print "B-R1: ".$b_tlxl->{R1_Qstart}."-".$b_tlxl->{R1_Qend}."\n" if defined $b_tlxl->{R1_Qstart};
    # print "B-R2: ".$b_tlxl->{R2_Qstart}."-".$b_tlxl->{R2_Qend}."\n" if defined $b_tlxl->{R2_Qstart};


    my $b_ref;
    sswitch ($tlx->{B_Rname}) {
      case "Breaksite": { $b_ref = $refs->{brk}; }
      case "Adapter": { $b_ref = $refs->{adpt}; }
      default: { $b_ref = $refs->{genome}; }
    }



    if ($i < $R2_idx || $R2_idx < 0) {
      # Segment earlier than the merge or no merge
      # Use R1
      $tlx->{B_Rstart} = $b_tlxl->{R1_Rstart};
      $tlx->{B_Rend} = $b_tlxl->{R1_Rend};
      $tlx->{B_Qstart} = $R1_map->[$b_tlxl->{R1_Qstart}-1];
      $tlx->{B_Qend} = $R1_map->[$b_tlxl->{R1_Qend}-1];

      $tlx->{B_R1_ID} = $b_tlxl->{R1_ID};

      my $Rseq = $b_ref->seq($tlx->{B_Rname},$tlx->{B_Rstart},$tlx->{B_Rend});
      $Rseq = reverseComplement($Rseq) unless $tlx->{B_Strand} == 1;
      # print "retrieved rseq ".$tlx->{B_Rname}." ".$tlx->{B_Rstart}." ".$tlx->{B_Rend}."\n$Rseq\n";

      $tlx->{B_CigarA} = soft_clip_cigar($b_tlxl->{R1_CigarA});

    } elsif ($i == $R2_idx) {
      # This is the segment we merged on
      if ($tlx->{B_Strand} == 1) {
        # "+" strand alignment
        $tlx->{B_Rstart} = $b_tlxl->{R1_Rstart};
        $tlx->{B_Rend} = $b_tlxl->{R2_Rend};

      } else {
        # "-" strand alignment
        $tlx->{B_Rstart} = $b_tlxl->{R2_Rstart};
        $tlx->{B_Rend} = $b_tlxl->{R1_Rend};
      }
      $tlx->{B_Qstart} = $R1_map->[$b_tlxl->{R1_Qstart}-1];
      $tlx->{B_Qend} = $R2_map->[$b_tlxl->{R2_Qend}-1];
      $tlx->{B_R1_ID} = $b_tlxl->{R1_ID};
      $tlx->{B_R2_ID} = $b_tlxl->{R2_ID};
      $tlx->{B_CigarA} = soft_clip_cigar($Cigar_array_ref);

    } else {
      # Segment post merge
      # Use R2
      $tlx->{B_Rstart} = $b_tlxl->{R2_Rstart};
      $tlx->{B_Rend} = $b_tlxl->{R2_Rend};
      $tlx->{B_Qstart} = $R2_map->[$b_tlxl->{R2_Qstart}-1];
      $tlx->{B_Qend} = $R2_map->[$b_tlxl->{R2_Qend}-1];
      $tlx->{B_R2_ID} = $b_tlxl->{R2_ID};

      my $Rseq = $b_ref->seq($tlx->{B_Rname},$tlx->{B_Rstart},$tlx->{B_Rend});
      $Rseq = reverseComplement($Rseq) unless $tlx->{B_Strand} == 1;
      # print "retrieved rseq ".$tlx->{B_Rname}." ".$tlx->{B_Rstart}." ".$tlx->{B_Rend}."\n$Rseq\n";

      $tlx->{B_CigarA} = soft_clip_cigar($b_tlxl->{R2_CigarA});

    }

    $tlx->{B_Cigar} = cigar_array_to_string($tlx->{B_CigarA});

    
    if ($i+1 <= $#$tlxls) {
      my $tlxl = $tlxls->[$i+1];

      # print "R1: ".$tlxl->{R1_Qstart}."-".$tlxl->{R1_Qend}."\n" if defined $tlxl->{R1_Qstart};
      # print "R2: ".$tlxl->{R2_Qstart}."-".$tlxl->{R2_Qend}."\n" if defined $tlxl->{R2_Qstart};

      $tlx->{Rname} = $tlxl->{Rname};
      $tlx->{Strand} = $tlxl->{Strand};


      my $ref;
      sswitch ($tlx->{Rname}) {
        case "Breaksite": { $ref = $refs->{brk}; }
        case "Adapter": { $ref = $refs->{adpt}; }
        default: { $ref = $refs->{genome}; }
      }


      if ($i+1 < $R2_idx || $R2_idx < 0) {
        $tlx->{Rstart} = $tlxl->{R1_Rstart};
        $tlx->{Rend} = $tlxl->{R1_Rend};
        $tlx->{Qstart} = $R1_map->[$tlxl->{R1_Qstart}-1];
        $tlx->{Qend} = $R1_map->[$tlxl->{R1_Qend}-1];

        my $Rseq = $ref->seq($tlx->{Rname},$tlx->{Rstart},$tlx->{Rend});
        $Rseq = reverseComplement($Rseq) unless $tlx->{Strand} == 1;
        # print "retrieved rseq ".$tlx->{Rname}." ".$tlx->{Rstart}." ".$tlx->{Rend}."\n$Rseq\n";

        $tlx->{CigarA} = soft_clip_cigar($tlxl->{R1_CigarA});

        $tlx->{R1_ID} = $tlxl->{R1_ID};

      } elsif ($i+1 == $R2_idx) {
        if ($tlx->{Strand} == 1) {
          $tlx->{Rstart} = $tlxl->{R1_Rstart};
          $tlx->{Rend} = $tlxl->{R2_Rend};
        } else{
          $tlx->{Rstart} = $tlxl->{R2_Rstart};
          $tlx->{Rend} = $tlxl->{R1_Rend};
        }
        $tlx->{Qstart} = $R1_map->[$tlxl->{R1_Qstart}-1];
        $tlx->{Qend} = $R2_map->[$tlxl->{R2_Qend}-1];

        $tlx->{CigarA} = soft_clip_cigar($Cigar_array_ref);

        $tlx->{R1_ID} = $tlxl->{R1_ID};
        $tlx->{R2_ID} = $tlxl->{R2_ID};

      } else {
        $tlx->{Rstart} = $tlxl->{R2_Rstart};
        $tlx->{Rend} = $tlxl->{R2_Rend};
        $tlx->{Qstart} = $R2_map->[$tlxl->{R2_Qstart}-1];
        $tlx->{Qend} = $R2_map->[$tlxl->{R2_Qend}-1];


        my $Rseq = $ref->seq($tlx->{Rname},$tlx->{Rstart},$tlx->{Rend});
        $Rseq = reverseComplement($Rseq) unless $tlx->{Strand} == 1;
        # print "retrieved rseq ".$tlx->{Rname}." ".$tlx->{Rstart}." ".$tlx->{Rend}."\n$Rseq\n";

        $tlx->{CigarA} = soft_clip_cigar($tlxl->{R2_CigarA});


        $tlx->{R2_ID} = $tlxl->{R2_ID};

      }

      $tlx->{Cigar} = cigar_array_to_string($tlx->{CigarA});

      if ($tlx->{Strand} == 1) {
        $tlx->{Junction} = $tlx->{Rstart};
      } else {
        $tlx->{Junction} = $tlx->{Rend};
      }
      
      $tlx->{J_Seq} = "";

    }

    push(@tlxs,$tlx);

    last if $i+1 >= $#$tlxls;
  }

  return(\@tlxs);

}

sub merge_alignments ($$) {

  my $tlxls = shift;
  my $refs = shift;

  my $R2_idx = firstidx {defined $_->{R2_Seq}} @$tlxls;

  if ($R2_idx < 0) {
    # Return a blank seq and some empty arraryrefs if not
    return("",[],[],[]);
  }

  # Pull the element to merge on out of the tlxl array
  my $tlxl = $tlxls->[$R2_idx];


  # We're gonna spend a while getting setup to do the merge
  # Here's the deal, in order to treat "+" and "-" strands the same
  # we have to take advantage of the symmetry
  # We do this by swapping reads 1 and 2 if its a "-" strand alignment
  # More work now and at the end when we have to flip "-" strands,
  # but way less work when we are actually ticking through and merging the reads

  my $Strand = $tlxl->{Strand};
  my $Rname = $tlxl->{Rname};
  my $Rstart1;
  my $Rend1;
  my $Rstart2;
  my $Rend2;

  my @Qseq1;
  my @Qseq2;
  my @Qual1;
  my @Qual2;
  my @Cigar1 = ();
  my @Cigar2 = ();
  my $Qstart1;
  my $Qend1;
  my $Qstart2;
  my $Qend2;

  # Need to know which sam object to use
  my $ref;
  sswitch ($Rname) {
    case "Breaksite": { $ref = $refs->{brk}; }
    case "Adapter": { $ref = $refs->{adpt}; }
    default: { $ref = $refs->{genome}; }
  }



  if ($Strand == 1) {
    # print "+ strand alignment\n";
    # Query sequence and quality info for each pair split into arrays
    @Qseq1 = split("",$tlxl->{R1_Seq});
    @Qseq2 = split("",$tlxl->{R2_Seq});
    @Qual1 = @{$tlxl->{R1_Qual}};
    @Qual2 = @{$tlxl->{R2_Qual}};

    # Create cigar arrays

    @Cigar1 = @{ expand_cigar_array(soft_clip_cigar($tlxl->{R1_CigarA})) };
    @Cigar2 = @{ expand_cigar_array(soft_clip_cigar($tlxl->{R2_CigarA})) };

    # foreach my $i (@{$tlxl->{R1_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar1, ($i->[0]) x $i->[1]);
    # }
    # foreach my $i (@{$tlxl->{R2_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar2, ($i->[0]) x $i->[1]);
    # }

    # Pretty straight forward this time through
    $Qstart1 =  $tlxl->{R1_Qstart};
    $Qend1 = $tlxl->{R1_Qend};
    $Qstart2 = $tlxl->{R2_Qstart};
    $Qend2 = $tlxl->{R2_Qend};

    $Rstart1 =  $tlxl->{R1_Rstart};
    $Rend1 = $tlxl->{R1_Rend};
    $Rstart2 = $tlxl->{R2_Rstart};
    $Rend2 = $tlxl->{R2_Rend};


  } else {
    # print "- strand alignment\n";
    # Query info again, this time flipping R1/R2, reversing (and complementing)
    @Qseq1 = split("",reverseComplement($tlxl->{R2_Seq}));
    @Qseq2 = split("",reverseComplement($tlxl->{R1_Seq}));
    @Qual1 = reverse @{$tlxl->{R2_Qual}};
    @Qual2 = reverse @{$tlxl->{R1_Qual}};

    # Create cigar arrays, same as before, but reversed

    @Cigar1 = reverse @{ expand_cigar_array(soft_clip_cigar($tlxl->{R2_CigarA})) };
    @Cigar2 = reverse @{ expand_cigar_array(soft_clip_cigar($tlxl->{R1_CigarA})) };

    # foreach my $i (reverse @{$tlxl->{R2_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar1, ($i->[0]) x $i->[1]);
    # }
    # foreach my $i (reverse @{$tlxl->{R1_CigarA}}) {
    #   next if $i->[0] eq "S";
    #   push(@Cigar2, ($i->[0]) x $i->[1]);
    # }

    # This time flip all R1 with R2
    # And since we're reverse complimenting
    $Qstart1 = $tlxl->{R2_Qlen} - $tlxl->{R2_Qend} + 1;
    $Qend1 = $tlxl->{R2_Qlen} - $tlxl->{R2_Qstart} + 1;
    $Qstart2 = $tlxl->{R1_Qlen} - $tlxl->{R1_Qend} + 1;
    $Qend2 = $tlxl->{R1_Qlen} - $tlxl->{R1_Qstart} + 1;

    $Rstart1 =  $tlxl->{R2_Rstart};
    $Rend1 = $tlxl->{R2_Rend};
    $Rstart2 = $tlxl->{R1_Rstart};
    $Rend2 = $tlxl->{R1_Rend};
  
  }

  # Retrieve reference sequence
  my @Rseq = split("",uc($ref->seq($Rname,$Rstart1,$Rend2)));


  # I think we're ready to merge!

  # Here's where we'll save the returned variables
  my @Qseq;
  my @Qmap1 = (0) x @Qseq1;
  my @Qmap2 = (0) x @Qseq2;
  my @Cigar;


  # We will be interating through these
  my $Qpos1 = $Qstart1;
  my $Qpos2 = $Qstart2;
  my $Rpos = $Rstart1;


  # Handle everything on R1 up to Qstart
  # Note!!! For "-" alignments, R1 is actually temporarily R2
  # (but we can forget that because symmetry)
  if ($Qstart1 > 1) {
    # print "premerge\n";

    my @R1_tail = 0..($Qstart1-2);
    push(@Qseq, @Qseq1[@R1_tail]);
    push(@Cigar,("S") x @R1_tail);
    @Qmap1[@R1_tail] = 1..@R1_tail;
  }


  # Slightly tricky... if there is a dovetail at the near end 
  # R1    |------------->
  # R2  <------------------------|
  # We have to push forward on the R2 query so it's caught up
  # Basically this means shifting off the front of the cigar array
  # and incrementing $Qpos2 until we're at the same Ref position as R1
  if ($Rstart2 < $Rstart1) {
    # print "pushing R2 premerge\n";
    my $Rpos_tmp = $Rstart2;
    while ($Rpos_tmp < $Rstart1) {
      my $c2 = shift(@Cigar2);
      sswitch ($c2) {
        case 'M':
        case 'X': {
          $Qpos2++;
          $Rpos_tmp++;
        }
        case 'D': {
          $Rpos_tmp++;
        }
        case 'I': {
          $Qpos2++;
        }
      }
    }
  }

  # Push ahead on R1 where there is no overlap
  # So either until the end of R1 or the start of R2
  # print "pushing R1 at $Rpos\n";
  while ($Rpos <= $Rend1 && $Rpos < $Rstart2) {
    my $c1 = shift(@Cigar1);
    sswitch ($c1) {
      case 'M':
      case 'X': {
        push(@Qseq,$Qseq1[$Qpos1-1]);
        if ($Qseq1[$Qpos1-1] eq $Rseq[$Rpos-$Rstart1]) {
          push(@Cigar,"M");
        } else {
          push(@Cigar,"X");
        }
        $Qmap1[$Qpos1-1] = scalar @Qseq;
        $Qpos1++;
        $Rpos++;
      }
      case 'D': {
        push(@Cigar,"D");
        $Rpos++;
      }
      case 'I': {
        push(@Qseq,$Qseq1[$Qpos1-1]);
        push(@Cigar,"I");
        $Qmap1[$Qpos1-1] = scalar @Qseq;
        $Qpos1++;
      }
    }
  }


  # Handles the gap situation
  # R1   |------------>
  # R2                      <--------------|
  if ($Rstart2 - $Rend1 > 1) {
    # print "handling gap at $Rpos\n";
    my $Gap_start = $Rend1 - $Rstart1 + 1;
    my $Gap_end = $Rstart2 - $Rstart1 - 1;
    push(@Qseq, map {lc $_} @Rseq[$Gap_start..$Gap_end]);
    push(@Cigar,("N") x ($Gap_end-$Gap_start+1));
    $Rpos = $Rstart2;
  }

  # Handles the overlap situation
  # R1  |------------->
  # R2       <--------------|
  if ($Rend1 >= $Rstart2) {
    # print "handling overlap at $Rpos\n";
    while ($Rpos <= $Rend1 && $Rpos <= $Rend2) {
      my $c1 = shift(@Cigar1);
      my $c2 = shift(@Cigar2);
      my $q1 = $Qseq1[$Qpos1-1];
      my $q2 = $Qseq2[$Qpos2-1];
      my $r = $Rseq[$Rpos-$Rstart1];


      sswitch ($c1) {
        case 'M':
        case 'X': {
          sswitch ($c2) {
            case 'M':
            case 'X': {
              # Match - Match
              # Take base if they match each other,
              # otherwise take the one that matches the reference
              # otherwise take the one with the higher qual score
              # Update maps and push forward on Q1 Q2 and R
              my $n;
              if ($q1 eq $q2) {
                $n = $q1;
              } elsif ($q1 eq $r || $q2 eq $r) {
                $n = $r;
              } else {
                $n = $Qual1[$Qpos1-1] > $Qual2[$Qpos2-1] ? $q1 : $q2;
              }
              push(@Qseq,$n);
              if ($n eq $r) {
                push(@Cigar,"M");
              } else {
                # print "$Rpos, $r, $Qpos1, $q1, $Qpos2, $q2\n";
                push(@Cigar,"X");
              }
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos1++;
              $Qpos2++;
              $Rpos++;
            }
            case 'D': {
              # Match - Del
              # Take matched base if it matches reference
              # Update maps, and push forward on Q1 and R
              if ($q1 eq $r) {
                push(@Qseq,$q1);
                push(@Cigar,"M");
              } else {
                push(@Cigar,"D");
              }
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qpos1++;
              $Rpos++;
            }
            case 'I': {
              # Match - Ins
              # Put cigar back for C1
              # Update maps and push forward on Q2
              unshift(@Cigar1,$c1);
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos2++;
            }
          }
        }
        case 'D': {
          sswitch ($c2) {
            case 'M':
            case 'X': {
              # Del - Match
              # Take matched base if it matches reference
              # Update maps, and push forward on Q2 and R
              if ($q2 eq $r) {
                push(@Qseq,$q2);
                push(@Cigar,"M");
              } else {
                push(@Cigar,"D");
              }
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos2++;
              $Rpos++;
            }
            case 'D': {
              # Del - Del
              # Push forward on R
              push(@Cigar,"D");
              $Rpos++;
            }
            case 'I': {
              # Del - Ins
              # Put cigar back for C1
              # Update maps and push forward on Q2
              unshift(@Cigar1,$c1);
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos2++;
            }
          }
        }
        case 'I': {
          sswitch ($c2) {
            case 'M':
            case 'X': {
              # Match - Ins
              # Put cigar back for C2
              # Update maps and push forward on Q1
              unshift(@Cigar2,$c2);
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qpos1++;
            }
            case 'D': {
              # Ins - Del
              # Put cigar back for C2
              # Update maps and push forward on Q1
              unshift(@Cigar2,$c2);
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qpos1++;
            }
            case 'I': {
              # Ins - Ins
              # Take higher quality base
              # Update maps and push forward on Q1 and Q2
              push(@Qseq, $Qual1[$Qpos1-1] > $Qual2[$Qpos2-1] ? $q1 : $q2);
              push(@Cigar,"I");
              $Qmap1[$Qpos1-1] = scalar @Qseq;
              $Qmap2[$Qpos2-1] = scalar @Qseq;
              $Qpos1++;
              $Qpos2++;
            }
          }
        }
      }
    }
  }


  # Push ahead on R2
  # print "pushing R2 at $Rpos\n";
  while ($Rpos <= $Rend2) {
    my $c2 = shift(@Cigar2);
    sswitch ($c2) {
      case 'M':
      case 'X': {
        push(@Qseq,$Qseq2[$Qpos2-1]);
        if ($Qseq2[$Qpos2-1] eq $Rseq[$Rpos-$Rstart1]) {
          push(@Cigar,"M");
        } else {
          push(@Cigar,"X");
        }
        $Qmap2[$Qpos2-1] = scalar @Qseq;
        $Qpos2++;
        $Rpos++;
      }
      case 'D': {
        push(@Cigar,"D");
        $Rpos++;
      }
      case 'I': {
        push(@Qseq,$Qseq2[$Qpos2-1]);
        push(@Cigar,"I");
        $Qmap2[$Qpos2-1] = scalar @Qseq;
        $Qpos2++;
      }
    }
  }

  # Same as before, if there's a dovetail at the far end
  # R1  |-------------------------->
  # R2       <-------------------|
  # Just clean up R1 query
  if ($Rend1 > $Rend2) {
    # I think actually we don't have to worry about this end
    # since we don't care about $Qpos1 and @Cigar1 by this point
  }

  # Handle everything on R2 after Qend
  if (scalar @Qseq2 > $Qend2) {
    # print "postmerge\n";
    my @R2_tail = ($Qend2)..$#Qseq2;
    push(@Qseq, @Qseq2[@R2_tail]);
    push(@Cigar,("S") x @R2_tail);
    @Qmap2[@R2_tail] = (scalar @Qseq - scalar @R2_tail - 1)..(scalar @Qseq);
  }

  my $Qseq = join("",@Qseq);


  # print "R1\n";
  # print join("",map {chr($_+33)} @Qual1)."\n";
  # print join("",@Qseq1)."\n";
  # print join(" ",@Qmap1)."\n";
  
  # print "R2\n";
  # print join("",map {chr($_+33)} @Qual2)."\n";
  # print join("",@Qseq2)."\n";
  # print join(" ",@Qmap2)."\n";  

  # print "merge\n";
  # print "$Qseq\n";


  # Flip things around if it was - strand alignment
  if ($Strand == -1) {
    $Qseq = reverseComplement($Qseq);
    @Cigar = reverse @Cigar;
    my @Qmap1_tmp = @Qmap1;
    my @Qmap2_tmp = @Qmap2;
    @Qmap1 = map {$_ == 0 ? 0 : (scalar @Qseq - $_ + 1) } reverse @Qmap2_tmp;
    @Qmap2 = map {$_ == 0 ? 0 : (scalar @Qseq - $_ + 1) } reverse @Qmap1_tmp;

    # print "reversed\n";
    # print join(" ",@Qmap1)."\n";
    # print join(" ",@Qmap2)."\n";  
    # print "$Qseq\n";
  }

  return($Qseq,\@Qmap1,\@Qmap2,compact_cigar_array(\@Cigar));

}


sub split_breaksite_reads ($) {
  my $read_obj = shift;
  my $params = $main::params;


  
  
}


sub find_random_barcode ($$) {
  my $read_obj = shift;
  my $barcode_length = shift;

  my $tlxs = $read_obj->{tlxs};
  my $R1_alns = $read_obj->{R1_alns};
  my $R2_alns = $read_obj->{R2_alns};


  my $barcode = "";


  if (defined $barcode_length && $barcode_length > 0) {
    
    debug_print("searching for barcode",2,$tlxs->[0]->{QnameShort});

    # Search through OCS first

    if (defined $tlxs->[$#$tlxs]->{Rname} && $tlxs->[$#$tlxs]->{Rname} eq "Adapter") {
      my $adapter_aln = $tlxs->[$#$tlxs];
      $barcode = substr($adapter_aln->{Seq},$adapter_aln->{Qstart} - $barcode_length - 1,$barcode_length);
    } else {
      
      my $adapter_aln;

      foreach my $R2_aln (values %$R2_alns) {
        if ($R2_aln->{Unmapped} == 0 && $R2_aln->{Rname} eq "Adapter" && $R2_aln->{Strand} == 1) {
          $adapter_aln = $R2_aln if ! defined $adapter_aln ||
                                    $R2_aln->{Qend} > $adapter_aln->{Qend};
        }
      }

      foreach my $R1_aln (values %$R1_alns) {
        if ($R1_aln->{Unmapped} == 0 && $R1_aln->{Rname} eq "Adapter" && $R1_aln->{Strand} == 1) {
          $adapter_aln = $R1_aln if ! defined $adapter_aln ||
                                    $R1_aln->{Qend} > $adapter_aln->{Qend};
        }
      }

      $barcode = substr($adapter_aln->{Seq},$adapter_aln->{Qstart} - $barcode_length - 1,$barcode_length) if defined $adapter_aln;

    }

  }

  foreach my $tlx (@$tlxs) {
    $tlx->{Barcode} = $barcode;
  }

}




1;
