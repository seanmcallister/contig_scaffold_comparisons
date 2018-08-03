#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::MoreUtils qw(uniq);

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Goals of script:
#Infiles: contigs.paths, scaffolds.paths, contigs.fasta, scaffolds.fasta, qc_contigs.fasta.
#To do:
#	1.	Match contig node names with edges --> contig dictionary.
#	2.	Match scaffold node names with edges --> scaffold dictionary.
#	3.	Identify scaffolds (at least 1 N).
#	4.	Change name of scaffold node to match the contig with the lowest node number.
#	5.	Print file showing scaffolds with >1 contig and their names.
#	6. 	Print file with old and new scaffold node names.


# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %options=();
getopts("c:f:s:p:q:m:l:h", \%options);

if ($options{h})
    {   print "\n\nHelp called:\nOptions:\n";
	print "-c = contigs.paths file from spades\n";
	print "-f = contigs.fasta file from spades\n";
	print "-q = QC'ed contig fasta file\n";
	print "-s = scaffolds.fasta file from spades\n";
	print "-p = scaffolds.paths file from spades\n";
	print "-m = sample header prefix for simplified headers\n";
	print "-l = minimum sequence length to export (Note: will export two files: one with all sequences and one to this minimum length)\n";
	print "-h = This help message\n\n";
	die;
    }

my %Contigs;
my %Scaffolds;
my %QCcontigs;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
&CONTIGread($options{f}, 1);
&SCAFFOLDread($options{s}, 2);
&QCCONTIGread($options{q}, 3);


################Open contig.paths file, match edge names with node names
###Contig.paths format =
#NODE_2_length_56311_cov_29.322
#36742+,64328+,54895-,66414-,66415+
#
#
###Corresponds to the following edges:
#EDGE_36742_length_#_cov_##
#EDGE_64328_length_#_cov_##
#EDGE_54895_length_#_cov_##'
#EDGE_66414_length_#_cov_##'
#EDGE_66415_length_#_cov_##
###NOTE: Length is always a whole number. Coverage is not always a decimal.
#
#
###Occasionally, get a semicolon separating paths with contigs that jump over a gap in the assembly graph:
#NODE_16_length_36540_cov_7.1333'
#23690-;
#30017-
#

#Open contigs.paths
print "   Reading contig paths file . . . \n";
open(IN2, "<$options{c}") or die "\n\nNADA $options{c} you FOOL!!!\n\n";
my @paths = <IN2>; close(IN2);
my $NODE = "NULL";
foreach my $line (@paths)
    {   chomp($line);
	if ($line =~ m/^NODE.+/)
		{	$NODE = $line;
			$NODE =~ s/\'$//;
		}
	if ($line =~ m/[+-]\;$/)
		{	my $regexline = $line;
			$regexline =~ s/(.+)\;$/$1/;
			my @edges;
			if ($regexline =~ m/,/) {@edges = split(',', $regexline);}
			else {push(@edges, $regexline);}
			foreach my $i (@edges)
				{	$i =~ s/([0-9]+)[+-]/$1/;
					$Contigs{$NODE}{'edges'} .= $i.",";	
				}	
		}
	if ($line =~ m/[0-9]+[+-]$/)
		{	my @edges;
			if ($line =~ m/,/) {@edges = split(',', $line);}
			else {push(@edges, $line);}
			foreach my $i (@edges)
				{	$i =~ s/([0-9]+)[+-]/$1/;
					$Contigs{$NODE}{'edges'} .= $i.",";
				}	
		}
    }

foreach my $i (sort keys %Contigs)
	{	my $edge_uniq = $Contigs{$i}{'edges'};
		my @edge_uniqs = split(',', $edge_uniq);
		my @uniq = uniq(@edge_uniqs);
		my @uniq_sorted = sort {$a <=> $b} @uniq;
		my $uniq = join(',', @uniq_sorted);
		$Contigs{$i}{'edges'} = $uniq;	
	}





#Open scaffolds.paths
print "   Reading scaffold paths file . . . \n";
open(IN3, "<$options{p}") or die "\n\nNADA $options{p} you FOOL!!!\n\n";
my @scaffpaths = <IN3>; close(IN3);
$NODE = "NULL";
foreach my $line (@scaffpaths)
    {   chomp($line);
	if ($line =~ m/^NODE.+/)
		{	$NODE = $line;
			$NODE =~ s/\'$//;
		}
	if ($line =~ m/[+-]\;$/)
		{	my $regexline = $line;
			$regexline =~ s/(.+)\;$/$1/;
			my @edges;
			if ($regexline =~ m/,/) {@edges = split(',', $regexline);}
			else {push(@edges, $regexline);}
			foreach my $i (@edges)
				{	$i =~ s/([0-9]+)[+-]/$1/;
					$Scaffolds{$NODE}{'edges'} .= $i.",";	
				}	
		}
	if ($line =~ m/[0-9]+[+-]$/)
		{	my @edges;
			if ($line =~ m/,/) {@edges = split(',', $line);}
			else {push(@edges, $line);}
			foreach my $i (@edges)
				{	$i =~ s/([0-9]+)[+-]/$1/;
					$Scaffolds{$NODE}{'edges'} .= $i.",";
				}		
		}
    }

foreach my $i (sort keys %Scaffolds)
	{	my $edge_uniq = $Scaffolds{$i}{'edges'};
		my @edge_uniqs = split(',', $edge_uniq);
		my @uniq = uniq(@edge_uniqs);
		my @uniq_sorted = sort {$a <=> $b} @uniq;
		my $uniq = join(',', @uniq_sorted);
		$Scaffolds{$i}{'edges'} = $uniq;	
	}

#Tag %Contigs with a QC_PASS = TRUE/FALSE
print "   Tagging contigs dictionary with QC pass/fail . . . \n";
foreach my $i (sort keys %Contigs)
	{	my $mod_contig_header = $Contigs{$i}{'HEAD'};
		$mod_contig_header =~ s/^(NODE_[0-9]+)_length.+/$1/;
		my $match_qc_header = $options{m}."_".$mod_contig_header;
		if (exists $QCcontigs{$match_qc_header})
			{	$Contigs{$i}{'QC_PASS'} = 'TRUE';
				#print "$i\tmatched\t$match_qc_header\n";
			}
		#else {print "$i failed QC\n";}
	}

#Identify scaffolds (at least 1 N). %Scaffolds MULTICONTIG = TRUE/FALSE.
print "   Identifying multicontig scaffolds . . . \n";
foreach my $i (sort keys %Scaffolds)
	{	my $seq = $Scaffolds{$i}{'gappy-ntseq'};
		if ($seq =~ m/N/)
			{	$Scaffolds{$i}{'MULTICONTIG'} = 'TRUE';}
	}

#Add new_header for scaffolds with only a single contig. Copy QC stats to these scaffolds.
print "   Adding new headers for single contig scaffolds . . . \n";
my @contigs_dict_keys;
my @contigs_dict_edges;
foreach my $j (sort keys %Contigs)
	{	push(@contigs_dict_keys, $j);
		push(@contigs_dict_edges, $Contigs{$j}{'edges'});
	}
foreach my $i (sort keys %Scaffolds)
	{  if ($Scaffolds{$i}{'MULTICONTIG'} eq 'FALSE')
		{   #print "Reading\t$i\n";
		    my $edge_test = $Scaffolds{$i}{'edges'};
		    foreach my $j (0..$#contigs_dict_edges)
			{	my $contigs_edge_test = $contigs_dict_edges[$j];
				if ($edge_test eq $contigs_edge_test)
					{	$Scaffolds{$i}{'new_header'} = $contigs_dict_keys[$j];
						$Scaffolds{$i}{'QC_PASS'} = $Contigs{$contigs_dict_keys[$j]}{'QC_PASS'};
						#print "$i\tmatches\t$Scaffolds{$i}{'new_header'}\n";
						#print "$edge_test\n";
						splice(@contigs_dict_keys, $j, 1);
						splice(@contigs_dict_edges, $j, 1);
						last;
					}
			}
		#    my @edge_matches = grep {$Contigs{$_}{'edges'} eq $edge_test} keys %Contigs;    
		#    if ($#edge_matches + 1 == 1)
		#	{	$Scaffolds{$i}{'new_header'} = $Contigs{$edge_matches[0]}{'HEAD'};
		#		$Scaffolds{$i}{'QC_PASS'} = $Contigs{$edge_matches[0]}{'QC_PASS'};
		#		print "$i\tmatches\t$Contigs{$edge_matches[0]}{'HEAD'}\n";
		#		print "$edge_test\n";
		#	}
		#    else {die "\n\nMORE THAN ONE EDGE MATCH:\n$i\n$edge_test\n<<@edge_matches>>\n\n";}
		}
	}

#Process multicontig scaffolds somehow. Interpret QC pass of any one contig in a scaffold = QC_PASS = TRUE for the scaffold.
print "   Processing multicontig scaffolds . . . \n";
foreach my $i (sort keys %Scaffolds)
	{	if ($Scaffolds{$i}{'MULTICONTIG'} eq 'TRUE')
			{	#print "Testing multicontig scaffold:  $i\n";
				my $seq = $Scaffolds{$i}{'gappy-ntseq'};
				my $AN = "AN";
				my @AN_count = $seq =~ /$AN/g;
				my $AN_count = $#AN_count + 1;
				#print "$AN_count\n";
				my $CN = "CN";
				my @CN_count = $seq =~ /$CN/g;
				my $CN_count = $#CN_count + 1;
				#print "$CN_count\n";
				my $TN = "TN";
				my @TN_count = $seq =~ /$TN/g;
				my $TN_count = $#TN_count + 1;
				#print "$TN_count\n";
				my $GN = "GN";
				my @GN_count = $seq =~ /$GN/g;
				my $GN_count = $#GN_count + 1;
				#print "$GN_count\n";
				my $NA = "NA";
				my @NA_count = $seq =~ /$NA/g;
				my $NA_count = $#NA_count + 1;
				#print "$NA_count\n";
				my $NC = "NC";
				my @NC_count = $seq =~ /$NC/g;
				my $NC_count = $#NC_count + 1;
				#print "$NC_count\n";
				my $NT = "NT";
				my @NT_count = $seq =~ /$NT/g;
				my $NT_count = $#NT_count + 1;
				#print "$NT_count\n";
				my $NG = "NG";
				my @NG_count = $seq =~ /$NG/g;
				my $NG_count = $#NG_count + 1;
				#print "$NG_count\n";
				my $total_count = $AN_count + $CN_count + $TN_count + $GN_count + $NA_count + $NC_count + $NT_count + $NG_count;
				my $total_NNNs = ($total_count)/2;
				my $total_contigs = $total_NNNs + 1;
				$Scaffolds{$i}{'total_contigs'} = $total_contigs;
				#print "Total contigs:  $Scaffolds{$i}{'total_contigs'}\n";
				my $scaf_edges = $Scaffolds{$i}{'edges'};
				my @edge_array = split(',', $scaf_edges);
				my @contig_matches;
				my @contig_qc_pass;
				foreach my $num (@edge_array)
					{	foreach my $d (sort keys %Contigs)
							{	my $contig_edges = $Contigs{$d}{'edges'};
								my @contig_edges = split(',', $contig_edges);
								foreach my $c_num (@contig_edges)
									{	if ($num == $c_num)
											{	push(@contig_matches, $Contigs{$d}{'HEAD'});
												push(@contig_qc_pass, $Contigs{$d}{'QC_PASS'});
											}
									}
							}
					}
				my @contig_call_uniq = uniq(@contig_matches);
				#print "Unique contig matches: @contig_call_uniq\n";
				#unless ($#contig_call_uniq + 1 == 1)
				if ($#contig_call_uniq + 1 != $Scaffolds{$i}{'total_contigs'})
					{	my $piss = $#contig_call_uniq + 1;
						print "WARNING: Scaffold does not match the predicted number of contigs: $i\tPredicted = $Scaffolds{$i}{'total_contigs'}\tActual = $piss\n";
						#DO SOMETHING DIFFERENT WITH THESE.
						#TRY TO MATCH THE CORRECT NUMBER OR SHOW USER AND HAVE THEM DECIDE!!
						my %Matched;
						foreach my $d (@contig_call_uniq)
							{	$Matched{$d}{'HEAD'} = $d;
								$Matched{$d}{'edges'} = $Contigs{$d}{'edges'};
								$Matched{$d}{'QC_PASS'} = $Contigs{$d}{'QC_PASS'};
								$Matched{$d}{'count'} = 0;
							}
						foreach my $d (@contig_matches)
							{	$Matched{$d}{'count'} += 1;}
						my @improper_match;
						foreach my $matched (sort keys %Matched)
							{	my $match_edges = $Matched{$matched}{'edges'};
								my @matched_edges = split(',', $match_edges);
								#print "<<<$match_edges>>>\n";
								my @match_truefalse;
								foreach my $hmmm (0..$#matched_edges)
									{ push(@match_truefalse, 'FALSE');
									}
								foreach my $s_edges (0..$#edge_array)
									{ foreach my $e_num (0..$#matched_edges)
										{ if ($edge_array[$s_edges] == $matched_edges[$e_num])
											{	$match_truefalse[$e_num] = 'TRUE';
											}
										}	
									}
								#print "<<<@match_truefalse>>>\n";
								foreach my $true_test (@match_truefalse)
									{	if ($true_test eq 'FALSE')
											{	push(@improper_match, $matched);	
											}	
									}
								#print "<<@improper_match>>\n";
							}
						my @uniq_improper_matches = uniq(@improper_match);	
						#print "<<<@uniq_improper_matches>>>\n";
						my @proper_matches;
						my @proper_truefalse;
						foreach (0..$#contig_call_uniq)
							{push(@proper_truefalse, 'TRUE');}
						foreach my $oranges (0..$#uniq_improper_matches)
							{	foreach my $apples (0..$#contig_call_uniq)
									{	if ($uniq_improper_matches[$oranges] eq $contig_call_uniq[$apples])
											{	$proper_truefalse[$apples] = 'FALSE';}
									}
							}
						foreach my $false_test (0..$#proper_truefalse)
							{	if ($proper_truefalse[$false_test] eq 'TRUE')
									{	push(@proper_matches, $contig_call_uniq[$false_test]);}
							}
						if ($#proper_matches + 1 == $Scaffolds{$i}{'total_contigs'})
							{	print "Scaffold to contig match fixed\n";
								print "Scaffold edges: $Scaffolds{$i}{'edges'}\n";
								print "Contig match: @proper_matches\n";
								print "Contig edges: ";
								my @check_multi_qc1;
								foreach my $stupid_check (@proper_matches)
									{	print "$Contigs{$stupid_check}{'edges'}\t";
										push(@check_multi_qc1, $Contigs{$stupid_check}{'QC_PASS'});
									}
								print "\n\n";
								foreach my $check_multi_qc1_pass (@check_multi_qc1)
									{	if ($check_multi_qc1_pass eq 'TRUE')
											{$Scaffolds{$i}{'QC_PASS'} = 'TRUE';}
									}
								my @sorted_multi_qc1_contig_matches = sort @proper_matches;
								$Scaffolds{$i}{'new_header'} = $sorted_multi_qc1_contig_matches[0];
								$Scaffolds{$i}{'multicontig_headers'} = '';		
								foreach my $g (sort @proper_matches)		
									{	 $Scaffolds{$i}{'multicontig_headers'} .= $g.",";}	
									
							}
						else { my $printactual = $#proper_matches + 1;
						       print "$i\tTHERE'S STILL A PROBLEM\tPredicted = $Scaffolds{$i}{'total_contigs'}\tActual = $printactual\n";
						       print "Let's get some human input here:\n";
						       print "Scaffold node edges: $Scaffolds{$i}{'edges'}\n";
						       print "Possible matches are numbered:\n";
						       foreach my $possible_matches (0..$#proper_matches)
								{	my $print_option = $possible_matches + 1;
									print "   Option $print_option: $proper_matches[$possible_matches]\t$Contigs{$proper_matches[$possible_matches]}{'edges'}\n";
								}
						       print "Input contig matches to keep - comma delimited (i.e. 1,4,5)\n";
						       my $user_choice = <STDIN>;
						       chomp($user_choice);
						       my @user_choice = split(',', $user_choice);
						       my @user_qc_match;
						       my @user_contig_match;
						       foreach my $user (@user_choice)
								{	my $hit = $user - 1;
									push(@user_contig_match, $proper_matches[$hit]);
									push(@user_qc_match, $Contigs{$proper_matches[$hit]}{'QC_PASS'});
								}
							foreach my $user_qc_pass (@user_qc_match)
								{  if ($user_qc_pass eq 'TRUE')
									{$Scaffolds{$i}{'QC_PASS'} = 'TRUE';}
								}
							my @sorted_user_contig_matches = sort @user_contig_match;
							$Scaffolds{$i}{'new_header'} = $sorted_user_contig_matches[0];
						        $Scaffolds{$i}{'multicontig_headers'} = '';
							foreach my $g (sort @user_contig_match)
							     {	 $Scaffolds{$i}{'multicontig_headers'} .= $g.",";}
							print "\n\n";
						     }
					
					
					}
				else
					{    #print "Successfull multicontig match\n";
					     foreach my $pass (@contig_qc_pass)
						     {	if ($pass eq 'TRUE')
								     {	$Scaffolds{$i}{'QC_PASS'} = 'TRUE';}
						     }
					     my @sorted_contig_matches = sort @contig_call_uniq;
					     $Scaffolds{$i}{'new_header'} = $sorted_contig_matches[0];
					     $Scaffolds{$i}{'multicontig_headers'} = '';
					     foreach my $g (sort @contig_call_uniq)
						     {	 $Scaffolds{$i}{'multicontig_headers'} .= $g.",";}
					}
				
			}
		
	}

#CHECK FOR NO NEW NAME AND DUPLICATES
print "   Checking for problems . . . \n";
my @test_dups;
foreach my $i (sort keys %Scaffolds)
	{	if ($Scaffolds{$i}{'new_header'} eq '')
			{	print "WARNING: $i does not have a new header match to a contig node name!\n";}
		unless ($Scaffolds{$i}{'new_header'} eq '')
			{	push(@test_dups, $Scaffolds{$i}{'new_header'});}
		#my @test_dups = grep {$Scaffolds{$_}{'new_header'} eq $Scaffolds{$i}{'new_header'}} keys %Scaffolds;
		#if (scalar @test_dups > 1)
		#	{	print "WARNING: Scaffold $i does not have a unique header, $Scaffolds{$i}{'new_header'}\n";}	
	}
my @uniq_test_dups = uniq(@test_dups);
if (scalar @test_dups != scalar @uniq_test_dups)
	{	print "WARNING: One or more scaffolds do not have unique new headers. Check the old/new node name log.\n";}


#Print file with scaffold names that were tossed due to QC.
print "   Saving outfiles . . . \n";
open(OUT1, ">scaffolds_not_passing_QC.txt");
print OUT1 "old_name\tnew_name\tmulticontig\tlength\n";
foreach my $s (sort keys %Scaffolds)
	{ if ($Scaffolds{$s}{'QC_PASS'} eq 'FALSE')
		{	print OUT1 "$Scaffolds{$s}{'HEAD'}\t$Scaffolds{$s}{'new_header'}\t$Scaffolds{$s}{'MULTICONTIG'}\t$Scaffolds{$s}{'N_removed_SIZE'}\n";	
		}	
	}
close(OUT1);

#Print file showing old and new scaffold names.
open(OUT2, ">scaffolds_old_new_names.txt");
print OUT2 "old_name\tnew_name\n";
foreach my $s (sort keys %Scaffolds)
	{	print OUT2 "$Scaffolds{$s}{'HEAD'}\t$Scaffolds{$s}{'new_header'}\n";}
close(OUT2);

#Print file with scaffolds >1 contig and the contigs that are in those scaffolds. And the QC stats of each contig and scaffold.
open(OUT3, ">scaffolds_multicontig_stats.txt");
print OUT3 "new_name\tscaffold_QC_status\tcontig_name\tcontig_QC_status\n";
foreach my $s (sort keys %Scaffolds)
	{	if ($Scaffolds{$s}{'MULTICONTIG'} eq 'TRUE')
		{	my $contig_list = $Scaffolds{$s}{'multicontig_headers'};
			my @contigs_list = split(',', $contig_list);
			my @qc_pass;
			foreach my $j (0..$#contigs_list)
				{	foreach my $i (sort keys %Contigs)
						{	if ($contigs_list[$j] eq $Contigs{$i}{'HEAD'})
								{push(@qc_pass, $Contigs{$i}{'QC_PASS'});}
						}
				}
			print OUT3 "$Scaffolds{$s}{'new_header'}\t$Scaffolds{$s}{'QC_PASS'}\t";
			foreach my $k (0..$#contigs_list)
				{	print OUT3 "$contigs_list[$k]\t$qc_pass[$k]\t";
				}
			print OUT3 "\n";
		}
	}
close(OUT3);

#Print fasta file with new scaffold node names: new and simplified.
#Print fasta file with new scaffold node names: QC scaffolds only, new and simplified.
#Print fasta file with new scaffold node names: QC scaffolds > min length only, new and simplified. 
open(OUT4, ">scaffolds_newnames.fasta");
open(OUT5, ">scaffolds_newnames_simplified.fasta");
open(OUT6, ">scaffolds_newnames_QCpass.fasta");
open(OUT7, ">scaffolds_newnames_simplified_QCpass.fasta");
open(OUT8, ">scaffolds_newnames_QCpass_atleast".$options{l}."bp.fasta");
open(OUT9, ">scaffolds_newnames_simplified_QCpass_atleast".$options{l}."bp.fasta");
foreach my $s (sort keys %Scaffolds)
	{	print OUT4 ">$Scaffolds{$s}{'new_header'}\n$Scaffolds{$s}{'gappy-ntseq'}\n";
		my $simplified_header = $Scaffolds{$s}{'new_header'};
		$simplified_header =~ s/^(NODE_[0-9]+)_length.+/$1/;
		print OUT5 ">".$options{m}."_".$simplified_header."\n$Scaffolds{$s}{'gappy-ntseq'}\n";
		if ($Scaffolds{$s}{'QC_PASS'} eq 'TRUE')
			{	print OUT6 ">$Scaffolds{$s}{'new_header'}\n$Scaffolds{$s}{'gappy-ntseq'}\n";
				print OUT7 ">".$options{m}."_".$simplified_header."\n$Scaffolds{$s}{'gappy-ntseq'}\n";
				if ($Scaffolds{$s}{'N_removed_SIZE'} >= $options{l})
					{	print OUT8 ">$Scaffolds{$s}{'new_header'}\n$Scaffolds{$s}{'gappy-ntseq'}\n";
						print OUT9 ">".$options{m}."_".$simplified_header."\n$Scaffolds{$s}{'gappy-ntseq'}\n";
					}
			}
	}
close(OUT4);
close(OUT5);
close(OUT6);
close(OUT7);
close(OUT8);
close(OUT9);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub CONTIGread
{	print "   Reading contig fasta file . . . \n";
	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
        my $filenumber = $_[1];
	open(IN, "<$infile") or die "\n\nNADA $infile 2 you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	#my $unid = $filenumber.1000000001;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		$Contigs{$data[0]}{'HEAD'}    = $data[0];       # store header
		$Contigs{$data[0]}{'QC_PASS'} = 'FALSE';
		$Contigs{$data[0]}{'edges'} = '';
		$Contigs{$data[0]}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		$Contigs{$data[0]}{'SIZE'}    = length($seq);   # store length
		$seq =~ s/\.//;
                $seq =~ s/\-//;
                $Contigs{$data[0]}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $Contigs{$data[0]}{'filenumber'} = $filenumber;
                #$unid += 1;
	}
	$/="\n";
}
sub SCAFFOLDread
{	print "   Reading scaffold fasta file . . . \n";
	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
        my $filenumber = $_[1];
	open(IN, "<$infile") or die "\n\nNADA $infile 2 you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	#my $unid = $filenumber.1000000001;                           # string to generate unique ids
	#open(OUTN, ">N_counts.txt");
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		$Scaffolds{$data[0]}{'HEAD'}    = $data[0];       # store header
		$Scaffolds{$data[0]}{'edges'} = '';
		$Scaffolds{$data[0]}{'new_header'} = '';
		$Scaffolds{$data[0]}{'MULTICONTIG'} = 'FALSE';
		$Scaffolds{$data[0]}{'QC_PASS'} = 'FALSE';
		$Scaffolds{$data[0]}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		$Scaffolds{$data[0]}{'SIZE'}    = length($seq);   # store length
		$seq =~ s/\.//;
                $seq =~ s/\-//;
                $Scaffolds{$data[0]}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $Scaffolds{$data[0]}{'filenumber'} = $filenumber;
                #my $N_count = $seq =~ tr/N/N/;
		#print OUTN "$N_count\n";
		$seq =~ s/N//;
		$Scaffolds{$data[0]}{'N_removed_SIZE'} = length($seq); 
		#$unid += 1;
	}
	#close(OUTN);
	$/="\n";
}
sub QCCONTIGread
{	print "   Reading qc fasta file . . . \n";
	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
        my $filenumber = $_[1];
	open(IN, "<$infile") or die "\n\nNADA $infile 2 you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	#my $unid = $filenumber.1000000001;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		$QCcontigs{$data[0]}{'HEAD'}    = $data[0];       # store header
		#$QCcontigs{$unid}{'edges'} = '';
		#$QCcontigs{$unid}{'gappy-ntseq'}   = uc($seq);       # store aligned sequence
		$QCcontigs{$data[0]}{'SIZE'}    = length($seq);   # store length
		$seq =~ s/\.//;
                $seq =~ s/\-//;
                #$QCcontigs{$unid}{'degapped-ntseq'} = uc($seq);     # store degapped sequence
                $QCcontigs{$data[0]}{'filenumber'} = $filenumber;
                #$unid += 1;
	}
	$/="\n";
}
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
