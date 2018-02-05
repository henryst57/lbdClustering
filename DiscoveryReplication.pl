#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

# usage: perl DiscoveryReplication.pl [vcluster_dir] [lbd_file] [vector_file] [cl_method] [target_cui]
# example: perl DiscoveryReplication.pl cluto-2.1.1/Linux/ data/rayFish_ltc_threshold_targetTermList data/1983_1985_window8 rb C0000545
#
# Input Parameters:
#  vcluster_dir - the directory containing the CLUTO vcluster executable
#  lbd_file - the file containing LBD output (target terms)
#  vector_file - the file containing word2vec vectors
#  cl_method - the clustering method to use: rb, rbr, direct, agglo, graph, 
#                                            bagglo
#
# Output:
#  Outputs to the (and if needed creates) the Results/ subdirectory
#    clustering info    - outputs the clusters file for each number of 
#                         clusters to: <target_cui>/<cl_method>/<numClusters> 
#                         the file contains the cluster number for the vector 
#                         at that line in the input vector file (also copied
#                         to the output directory as: <vector_file>_v). This
#                         file may be interpreted as vector at line X is in 
#                         cluster number at line X of the cluster file 
#    execution time     - appends the execution time to 
#                         ExecutionTimes/<target_cui>
#    target CUI Ranking - outputs the ranking of the <target_cui> at each 
#                         interval to: Rankings/<target_cui>_<cl_method>
#    visualization data - outputs files for each number of clusters that 
#                         contain the data needed by visualization program. 
#                         This file contains:
#                         CUI \t PreferredTerm \t ClusterName 
#                         and is output to:
#                      VisualizationData/<target_cui>/<cl_method>/<numClusters>
#
#

Main();


#########################################
# Main execution loop
sub Main {
    # initialize variables and get arguments
    my $start = time;
    my ($vcluster_dir, $lbd_file, $vector_file, $cl_method, $target_cui) = GetArgs();
    my $base_dir = cwd();	

    # set up output directories
    my $cl_soln_dir = CreateSolnDir($base_dir, $cl_method, $target_cui);
    my $ranking_dir = CreateRankingDir();
    my $vis_dir = CreateVisDir($target_cui, $cl_method);
    #my $tree_dir = CreateTreeDir();

    # read in data from file
    my ($og_cui_list, $cui_scores, $cui_terms) = ReadLBDData($lbd_file);
    my ($vectors, $new_cui_list) = ExtractVectors($og_cui_list, $vector_file);
    my $matrix_size = GetMatrixSize($vectors);

    # print updated data to file	
    my $v_ifile = PrintVClusterInputFile($vectors, $vector_file, $cl_soln_dir, $matrix_size);

    # run clustering	
    my $clustered_files = RunVCluster($base_dir, $cl_soln_dir, $vcluster_dir, $v_ifile, $cl_method, $matrix_size);

    # get rankings
    my %rankings;
    my %visualization_output;
    foreach my $num_clusters (sort {$a <=> $b} keys %$clustered_files){
	my $clustered_file = $clustered_files->{$num_clusters};
	my $clusters = ExtractClusters($clustered_file, $cl_soln_dir, $new_cui_list, $target_cui);

        # labelling of clusters for visualization
	my $centroids = CalculateCentroids($vectors, $clusters, $new_cui_list, $matrix_size);
	my $cluster_names = LabelClusters($centroids, $clusters, $vectors, $cui_terms, $new_cui_list);
	$visualization_output{$num_clusters} = CreateVisualizationOutput($cluster_names, $cui_terms, $clusters);

	my $target_cluster_i = SaveTargetCuiClusterIndex($clusters, $target_cui);
	my $cluster_scores = CalculateClusterScores($cui_scores, $clusters, $target_cluster_i);
	my $target_cui_intracluster_rank = RankTargetCuiWithinCluster($target_cluster_i, $target_cui, $clusters, $cui_scores);
	$rankings{$num_clusters} = RankTargetCui($cluster_scores, $target_cluster_i, $clusters, $target_cui_intracluster_rank);
    }

    # print rankings to file
    PrintClusteredRankings(\%rankings, $target_cui, $cl_method, $ranking_dir);
    PrintVisualizationOutput(\%visualization_output, $vis_dir);

    # get and print execution time
    my $execution_time = time - $start;
    PrintExecutionTime($execution_time, $target_cui, $cl_method);
}

sub PrintVisualizationOutput {
    my ($visualization_output, $vis_dir) = (@_);

    opendir my $dh, $vis_dir or die "Can't open $vis_dir: $!";
    foreach my $num_clusters (keys %$visualization_output){
	my $f = $vis_dir . '/' . $num_clusters;
	open my $fh, '>', $f or die "Can't open $f: $!";
	foreach my $line (@{$visualization_output->{$num_clusters}}){
	    print $fh "$line\n";
	}
	close $fh;
    }
    closedir $dh;
}

sub CreateVisualizationOutput {
    my ($cluster_names, $cui_terms, $clusters) = (@_);
    my @visualization_output;
    my $cluster_index = 0;
    foreach my $cluster_name (@$cluster_names){
	foreach my $cui (@{$clusters->{$cluster_index}}){
	    my $term = $cui_terms->{$cui};
	    my $str = "$cui\t$term\t$cluster_name";
	    push @visualization_output, $str;
	}
	$cluster_index += 1;
    }
    return \@visualization_output;
}

sub LabelClusters {
    my ($centroids, $clusters, $vectors, $cui_terms, $cui_list) = (@_);
    # labels cluster centroid with vector closest to the centroid, determined with cosine similarity 
    my @cluster_names;
    my %cui_list = reverse %$cui_list;
    foreach my $cluster_index (keys %$centroids){
	my $centroid_vector_index;
	my $max_cos_sim = -1; # start with minimum value for cosine similarity
	my @centroid = @{$centroids->{$cluster_index}};
	foreach my $cui (@{$clusters->{$cluster_index}}){
	    my $vector_index = $cui_list{$cui};
	    my @vector = @{$vectors->{$vector_index}};
	    my $cos_sim = CosSim(\@centroid, \@vector);
	    if ($cos_sim > $max_cos_sim){
		$max_cos_sim = $cos_sim; # update maximum cosine similarity value
		$centroid_vector_index = $vector_index;	# update centroid vector index
	    }
	}
	my $cui = $cui_list->{$centroid_vector_index};
	my $term = $cui_terms->{$cui};
	push @cluster_names, $term;
    }
    return \@cluster_names;
}

sub CosSim {
    # calculates cosine similarity between two vectors (aka the dot product)
    # cosine similarity forumla (Vector a, Vector b):
    # (a * b) / (|a|*|b|)
    my ($a, $b) = (@_);
    my @a = @$a;
    my @b = @$b;
    my $numerator;
    my $a_sq_sum;
    my $b_sq_sum;

    foreach my $i (0 .. $#a){
	$numerator += $a[$i]*$b[$i];
	$a_sq_sum += $a**2;
	$b_sq_sum += $b**2;
    }
    my $denominator = sqrt($a_sq_sum)*sqrt($b_sq_sum);
    my $cos_sim = $numerator/$denominator;
    return $cos_sim;
}

sub CalculateCentroids {
    # calculates the "centroid" (average value of all vectors) for all clusters for a particular clustering solution
    my ($vectors, $clusters, $cui_list, $matrix_size) = (@_);
    my %centroids;
    my %cui_list = reverse %$cui_list;
    my $num_columns = @$matrix_size[1];

    foreach my $cluster_index (keys %$clusters){ # iterate through clusters
	my $vectors_per_cluster = 0; # counter for vectors per cluster
	$centroids{$cluster_index} = [(0) x $num_columns]; # centroid values: initialize array of size num_columns to 0
	foreach my $cui (@{$clusters->{$cluster_index}}){ # for each cui within the cluster
	    $vectors_per_cluster++; # increment vectors per cluster
	    my $index = $cui_list{$cui};
	    my @vector = @{$vectors->{$index}}; # get vector values from vectors hash
	    foreach my $col (0 .. $num_columns-1){ # add vector to centroid values for cluster
		$centroids{$cluster_index}[$col] += $vector[$col]; # sum the values at each column
	    }
	}
	foreach my $col (0.. $num_columns-1){
	    $centroids{$cluster_index}[$col] /= $vectors_per_cluster; # take the average value over num vectors per cluster
	}
    }
    return \%centroids;
}

sub PrintExecutionTime {
    my ($execution_time, $target_cui, $cl_method) = (@_);
    my $execution_times_dir = 'Results/' . "ExecutionTimes";
    CreateDir($execution_times_dir);
    my $execution_times_file = $execution_times_dir . '/' . $target_cui;
    $execution_time /= 60;
    open my $fh, '>>', $execution_times_file or die "Can't open $execution_times_file: $!";
    print $fh "$cl_method\t$execution_time\n";
    close $fh;
}

sub PrintClusteredRankings {
    my ($rankings, $target_cui, $cl_method, $ranking_dir) = (@_);

    my $results_file = $ranking_dir . '/' . $target_cui . '_' . $cl_method;
    open my $fh, '>', $results_file, or die "Can't open $results_file: $!";
    foreach my $num_clusters (sort {$a <=> $b} keys %$rankings){
	my $rank = $rankings->{$num_clusters};
	print $fh "$rank\n";
    }
    close $fh;
}

sub ReverseClusterScores {
    my ($cluster_scores) = (@_);
    my %rev_cluster_scores;
    foreach my $cluster_i (keys %$cluster_scores){
	my $score = $cluster_scores->{$cluster_i};
	push @{$rev_cluster_scores{$score}}, $cluster_i;

    }
    return \%rev_cluster_scores;
}

sub RankTargetCui {
    # ranks the target cui
    # puts target cui rank as lowest cluster rank possible for all clusters of the same rank
    my ($cluster_scores, $target_cluster_i, $clusters, $target_cui_intracluster_rank) = (@_);
    my $rev_cluster_scores = &ReverseClusterScores($cluster_scores);
    my $target_cui_rank = 0;
    foreach my $score (sort {$b <=> $a} keys %$rev_cluster_scores){ # sort scores in descending order
	if (grep /$target_cluster_i/, @{$rev_cluster_scores->{$score}}){ # check if one of the clusters at that score contains target cui
	    foreach my $cl_i (@{$rev_cluster_scores->{$score}}){ # if it does, for each cluster,
		if ($cl_i != $target_cluster_i){ # that is not the target cluster,
		    my $num_cuis = @{$clusters->{$cl_i}}; # get the number of cuis in the cluster
		    $target_cui_rank += $num_cuis; # add this to the current target cui rank					
		} # this is how we rank the target cui cluster as the lowest ranked cluster for clusters of the same score
		else {
		    $target_cui_rank += $target_cui_intracluster_rank; # then, for the cluster that the target cui belongs to, we add its intracluster rank
		    return $target_cui_rank; # stop here and return the target cui
		}
	    }
	}
	else { # if cluster doesn't contain the target cui,
	    foreach my $cl_i (@{$rev_cluster_scores->{$score}}){
		my $num_cuis = @{$clusters->{$cl_i}}; # add up the number of cuis for all clusters at that score
		$target_cui_rank += $num_cuis; # and add this to the current target cui rank
	    }
	}			
    }
}

sub RankTargetCuiWithinCluster {
    my ($target_cluster_i, $target_cui, $clusters, $cui_scores) = (@_);
    my $target_cui_intracluster_rank;
    my $target_cui_score = $cui_scores->{$target_cui};
    my %temp;
    foreach my $cui (@{$clusters->{$target_cluster_i}}){
	my $score = $cui_scores->{$cui};
	push @{$temp{$score}}, $cui;
    }
    foreach my $score (sort {$b <=> $a} keys %temp){
	$target_cui_intracluster_rank += @{$temp{$score}};
	if ($score == $target_cui_score){
	    return $target_cui_intracluster_rank;
	}
    }
}

sub CalculateClusterScores {
    my ($cui_scores, $clusters) = (@_);	
    my %cluster_scores;

    my $cluster_index = 0;
    foreach my $cluster_index (keys %$clusters){
	$cluster_scores{$cluster_index} = 0;
	foreach my $cui (@{$clusters->{$cluster_index}}){
	    my $score = $cui_scores->{$cui};
	    $cluster_scores{$cluster_index} += $score;
	}
	$cluster_index++;
    }
    return \%cluster_scores;
}

sub SaveTargetCuiClusterIndex {
    my ($clusters, $target_cui) = (@_);
    foreach my $cluster_i (keys %$clusters){
	foreach my $cui (@{$clusters->{$cluster_i}}){
	    if ($cui eq $target_cui){
		return $cluster_i;
	    }
	}
    }
}

sub ExtractClusters {
    my ($clustered_file, $cl_soln_dir, $new_cui_list, $target_cui) = (@_);
    my %clusters;

    if (! -e "$cl_soln_dir/$clustered_file" && ! -f "$cl_soln_dir/$clustered_file") { return;}
    open my $fh, '<', "$cl_soln_dir/$clustered_file" or die "Can't open $cl_soln_dir/$clustered_file: $!"; # open file in read only mode
    my $index = 0; # cui index corresponds with line number
    while (my $line = <$fh>){ # read each line
	my $cluster_index = $line; # cluster index corresponds with line contents
	$cluster_index =~ s/\s+//; # clear excess whitespace
	my $cui = $new_cui_list->{$index}; # get cui
	push @{$clusters{$cluster_index}}, $cui; # array of clusters, index is cluster index, content is array of cuis
	$index++;
    }	
    close $fh;
    return \%clusters;
}

sub TargetCuiRankForOneCluster {
    my ($target_cui, $cui_scores) = (@_);
    my %temp_scores;
    my $target_cui_score;
    my $target_cui_rank;

    foreach my $cui (keys %$cui_scores){ # iterates through cui scores
	my $score = $cui_scores->{$cui};
	push @{$temp_scores{$score}}, $cui; # creates a temp_scores hash of arrays with score as key and array of cui's with this score
	if ($cui eq $target_cui){ # if cui is the target cui, saves the target cui score
	    $target_cui_score = $score;
	}
    }

    foreach my $score (sort {$b <=> $a} keys %temp_scores){ # iterates through temp_scores scores HOA in descending order
	$target_cui_rank += @{$temp_scores{$score}}; # adds number of cuis that have that score to the target cui rank
	if ($score == $target_cui_score){ # if score is the target cui's score,
	    return $target_cui_rank; # return the target cui rank.
	} # thus the target cui is ranked by score, and if it has the same score as other cui's,
    } # it is ranked lowest out of the group.
}

sub CreateDir {
    # checks if dir exists. if not, creates dir.
    my ($dir) = (@_);
    if (! -e $dir || ! -d $dir){
	mkdir $dir;
    }
}

sub RunVCluster {
    my ($base_dir, $cl_soln_dir, $vcluster_dir, $v_ifile, $cl_method, $matrix_size) = (@_);

    my %clustered_files;

    my $num_vectors = @$matrix_size[0];
    my $file_to_cluster = $base_dir . '/' . $cl_soln_dir . '/' . $v_ifile;
    for (my $i = 0; $i < 1.05; $i += 0.05){
	my $num_clusters = int($i * $num_vectors);

	print STDERR "i,num_clusters = $i, $num_clusters\n";
	
	if ($num_clusters == 0){ $num_clusters = 1; }
	my $clustered_filename = "$v_ifile.clustering.$num_clusters";
	$clustered_files{$num_clusters} = $clustered_filename;
	my $fp = $base_dir . '/' . $cl_soln_dir . '/' . $clustered_filename;
	my $tree_filename = $base_dir . '/' . $cl_soln_dir . '/' . $clustered_filename . '_tree';
	if (! -e $fp){
	    #my $cmd = "./vcluster -clmethod=$cl_method $file_to_cluster -showtree -cltreefile $tree_filename $num_clusters";
	     my $cmd = "./vcluster -clmethod=$cl_method $file_to_cluster -showtree $num_clusters";
	    print "$cmd\n";
	    chdir $vcluster_dir;
	    my $cluster_cmd = `$cmd`;
	    chdir $base_dir;
	} else {
	    die ("Error: clustering output already exists: $fp\n");
	}
    }
    return \%clustered_files;
}

sub PrintVClusterInputFile {
    # print new vector file with _v appended
    my ($vectors, $vector_file, $cl_soln_dir, $matrix_size) = (@_);

    my $num_rows = @$matrix_size[0];
    my $num_columns = @$matrix_size[1];

    $vector_file  =~ s/.*\/(.*?)/$1/;
    my $v_ifile = $vector_file . '_v';
    my $v_ifile_path = $cl_soln_dir . '/' . $v_ifile;

    if (-e $v_ifile_path && -f $v_ifile_path){
	return $v_ifile;
    }

    open my $fh, '>', "$v_ifile_path" or die "Can't open $v_ifile_path: $!"; # open file in create/write/truncate mode
    print $fh "$num_rows $num_columns\n";
    foreach my $vector_i (sort {$a <=> $b} keys %$vectors){
	print $fh join (' ', @{$vectors->{$vector_i}});
	print $fh "\n";
    }
    close $fh;

    return $v_ifile;
}

sub GetMatrixSize {
    my ($vectors) = (@_);
    my @matrix_size;
    $matrix_size[0] = keys %$vectors;	# number of rows
    $matrix_size[1] = @{$vectors->{0}};	# number of columns
    return \@matrix_size;
}

sub ExtractVectors {
    my ($og_cui_list, $vector_file) = (@_);
    my %vectors; # later, sort vectors by descending rank and print
    my %new_cui_list;
    my $index = 0;
    open my $fh, '<', "$vector_file" or die "Can't open $vector_file: $!";
    while (my $line = <$fh>){
	if ($line =~ /^(C\d{7})(.+)/){ # if line contains a cui/vector pair
	    my $cui = $1;		
	    if (exists $og_cui_list->{$cui}){ # if exists in target terms list (not all target terms will be in the vector list due to w2v threshold)
		my $vector = $2;
		my @vector_vals = ($vector =~ /(-?\d.\d+)/g);
		$vectors{$index} = [@vector_vals];
		$new_cui_list{$index} = $cui;
		$index += 1;
	    }
	}
    }
    close $fh;
    return \%vectors, \%new_cui_list;
}

sub ReadLBDData {
    my ($lbd_file) = (@_);
    my %og_cui_list;
    my %cui_scores;
    my %cui_terms;
    open my $fh, '<', "$lbd_file" or die "Can't open $lbd_file: $!";
    while (my $line = <$fh>) {
	if ($line =~ /(\d+)\t(\d+.\d+)\t(C\d{7})\t(.+)/){
	    my $index = $1 - 1;
	    my $score = $2;
	    my $cui = $3;
	    my $term = $4;
	    $og_cui_list{$cui} = $index;
	    $cui_scores{$cui} = $score;
	    $cui_terms{$cui} = $term;
	}
    }
    close $fh;
    if (keys %og_cui_list == 0){
	print "Invalid data in ltc file: $lbd_file.\n";
	exit;
    }
    return \%og_cui_list, \%cui_scores, \%cui_terms;
}

sub CreateVisDir {
    my ($target_cui, $cl_method) = (@_);

    my $vis_dir = 'Results/VisualizationData';
    CreateDir($vis_dir);

    my $target_cui_vis_dir = $vis_dir . '/' . $target_cui;
    CreateDir($target_cui_vis_dir);

    my $cl_method_vis_dir = $target_cui_vis_dir . '/' . $cl_method;
    CreateDir($cl_method_vis_dir);

    print "Visualization data will be printed to $cl_method_vis_dir\n";
    return $cl_method_vis_dir;
}

sub CreateRankingDir {
    my $ranking_dir = 'Results/Rankings';
    CreateDir($ranking_dir);
    print "Target term rankings will be printed to $ranking_dir\n";
    return $ranking_dir;
}

sub CreateSolnDir {
    my ($base_dir, $cl_method, $target_cui) = (@_);

    my $results_dir = "Results";
    my $target_cui_dir = $results_dir .'/' . $target_cui;
    my $cl_soln_dir = $target_cui_dir . '/' . $cl_method;

    &CreateDir($results_dir);
    &CreateDir($target_cui_dir);
    &CreateDir($cl_soln_dir);

    print "Clustering solutions will be printed to: $cl_soln_dir\n";
    return $cl_soln_dir;
}

sub GetArgs {
    my $vcluster_dir = $ARGV[0];
    my $lbd_file = $ARGV[1];
    my $vector_file = $ARGV[2];
    my $cl_method = $ARGV[3];
    my $target_cui = $ARGV[4];
    &CheckArgv($vcluster_dir, $lbd_file, $vector_file, $cl_method, $target_cui);
    return $vcluster_dir, $lbd_file, $vector_file, $cl_method, $target_cui;
}

sub CheckArgv {
    my ($vcluster_dir, $lbd_file, $vector_file, $cl_method, $target_cui) = (@_);
    &CheckNumArgs();
    &CheckVcluster($vcluster_dir);
    &CheckCLMethod($cl_method);
    &CheckFileErr($lbd_file);
    &CheckFileErr($vector_file);
    &CheckCui($target_cui);
}

sub CheckNumArgs {
    if ($#ARGV < 4) {
	print "Incorrect number of arguments. Program requires 5 arguments to run.\n";
	print "Usage: perl DiscoveryReplication.pl [vcluster_dir] [ltc_file] [vector_file] [cl_method] [target_cui]\n";
	exit;
    }
}

sub CheckVcluster {
    my ($vcluster_dir) = (@_);
    my $vcluster = $vcluster_dir . '/vcluster';
    if (! -e "$vcluster_dir/vcluster" || ! -d $vcluster_dir || ! -e $vcluster){
	print "Incorrect directory for vcluster program: $vcluster_dir.\n";
	exit;
    }
}

sub CheckCLMethod {
    my ($cl_method) = (@_);
    my @valid_cl_methods = ('rb', 'rbr', 'direct', 'agglo', 'graph', 'bagglo');
    if (!(grep /$cl_method/, @valid_cl_methods)) {
	print "Invalid method. Supported options are 'rb', 'rbr', 'direct', 'agglo', 'graph', 'bagglo'.\n";
	print "See CLUTO manual for descriptions of each clustering method (-clmethod).\n";
	exit;
    }
}

sub CheckFileErr {
    my ($file) = (@_);
    if (! -e $file || ! -f $file){
	print "File not found: $file.\n";
	exit;
    }
}

sub CheckCui {
    my ($cui) = (@_);
    if (!($cui =~ /^C\d{7}$/)){
	print "Incorrect cui format: $cui. Correct format is C followed by 7 digits, eg. C0000545.\n";
	exit;
    }
}
