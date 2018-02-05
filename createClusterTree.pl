#!/usr/bin/perl
use strict;
use warnings;

# usage: perl createClusterTree.pl [vcluster_dir] [lbd_file] [vector_file] [cl_method] [out_file]
# example: perl createClusterTree.pl ../cluto-2.1.1/Linux/ ../data/rayFish_ltc_threshold_targetTermList ../data/1983_1985_window8 rb output/rayFish_clusterTree
#
# Input Parameters:
#  vcluster_dir - the directory containing the CLUTO vcluster executable
#  lbd_file - the file containing LBD output (target terms)
#  vector_file - the file containing word2vec vectors
#  cl_method - the clustering method to use: rb, rbr, direct, agglo, graph, 
#                                            bagglo
#  out_file - the file to output the cluster tree to
#
#
# Output:
#  
#
#



Main();


########################################
# Begin Code
########################################

#main routine to perform clustering and output
sub Main {

    #get the input parameters
    my $output_dir = 'output/';
    &CreateDir($output_dir);
    my ($vcluster_dir, $lbd_file, $vector_file, $cl_method, $out_file) = GetArgs();
   
    # read in data from file
    print STDERR "   Reading Input File\n";
    my ($og_cui_list, $cui_scores, $cui_terms) = ReadLBDData($lbd_file);
    my ($vectors, $cui_list) = ExtractVectors($og_cui_list, $vector_file);
    my $matrix_size = GetMatrixSize($vectors);

    # run clustering	
    print STDERR "   Clustering\n";
    my $v_ifile = PrintVClusterInputFile($vectors, $vector_file, $output_dir, $matrix_size);
    (my $clusters_file, my $tree_file) = RunVCluster(
	$output_dir, $vcluster_dir, $v_ifile, $cl_method, $matrix_size);

    #cluster labelling and scoring
    print STDERR "Labelling Clusters\n";
    my $clusters = ExtractTree($tree_file, $cui_list);
    my $centroids = CalculateCentroids($vectors, $clusters, $cui_list, $matrix_size);
    my $cluster_names = LabelClusters($centroids, $clusters, $vectors, $cui_terms, $cui_list);

    #rank the clusters
    print STDERR "Calcualting Cluster Scores\n";
    my $cluster_scores = CalculateClusterScores($cui_scores, $clusters);

    #print out the cluster info
    print STDERR "Outputting Results\n";
    open OUT, ">$out_file" or die ("ERROR: cannot open out_file: $out_file\n");
    foreach my $clusterID (keys %{$clusters}) {
	print OUT "$clusterID - ${$cluster_names}{$clusterID} - ${$cluster_scores}{$clusterID}: ";
	foreach my $cui (@{${$clusters}{$clusterID}}) {
	    print OUT "$cui, "
	}
	print OUT "\n";
    }

    print STDERR "Done\n";
}



########################################
# Cluster Scoring
########################################
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



########################################
# Cluster Extracting
########################################
sub ExtractTree {
    my $tree_file = shift;
    my $cui_list = shift;
    
    #read the parent tree file
    open IN, $tree_file or die ("ERROR: cannot open tree_file: $tree_file\n");
    my %parents = ();
    my $index = 0;
    while (my $line = <IN>) {
	my ($clusterNum, $val1, $val2) = split(/ /, $line);
	$parents{$index} = $clusterNum;
	$index++;
    }
    close IN;

    #convert the parent tree into a children list
    #initialize the lsit
    my %children = ();
    foreach my $id (keys %parents) {
	my @emptyArray = ();
	$children{$id} = \@emptyArray;
    }
    #find the children of each id (cluster number ID)
    foreach my $id1 (keys %parents) {
	#find all children of $id1
	foreach my $id2 (keys %parents) {
	    #see if the parent of $id2 is $id1
	    if ($parents{$id2} == $id1) {
		#add $id2 to the list of children for $id1
		push @{$children{$id1}}, $id2;
	    }
	}
    }

    #convert the children list to a list of concepts in each cluster
    my %clusters = ();
    foreach my $id (keys %children) {
	$clusters{$id} = &getChildCuiIDs(\%children, $id);
    }

    #convert the CUI IDs to CUIs
    foreach my $id (keys %clusters) {
	my @cuis = ();
	foreach my $cuiID (@{$clusters{$id}}) {
	    push @cuis, ${$cui_list}{$cuiID};
	}
	$clusters{$id} = \@cuis;
    }

    #return the hash of clusters and their CUIs
    return \%clusters
}


#Recursively gets the descendent leaf nodes (CUIs)
# if the children tree. Returns an array ref of CUI IDs
sub getChildCuiIDs {
    my $childrenRef = shift;
    my $id = shift;

    #intialize
    my @children = @{${$childrenRef}{$id}}; 
    my @CUIs = ();

    #base case, this is a leaf node (no children)
    if (scalar @children == 0) {
	push @CUIs, $id;
	return \@CUIs;
    }

    #else keep going down until a leaf node is hit
    foreach my $childID(@children) {
	push @CUIs, @{&getChildCuiIDs($childrenRef, $childID)};
    }
    return \@CUIs;
}



########################################
# Cluster Labelling
########################################

#calculates the "centroid" (average value of all vectors) for all 
# clusters for a particular clustering solution
sub CalculateCentroids {
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

#labels cluster centroid with vector closest to the centroid, 
# determined with cosine similarity 
sub LabelClusters {
    my ($centroids, $clusters, $vectors, $cui_terms, $cui_list) = (@_);
   
    my %cluster_names = ();;
    my %cui_list = reverse %$cui_list;
    foreach my $cluster_index (keys %$centroids){
	my $centroid_vector_index;
	my $max_cos_sim = -1; # start with minimum value for cosine similarity
	my @centroid = @{$centroids->{$cluster_index}};
	foreach my $cui (@{$clusters->{$cluster_index}}){
	    my $vector_index = $cui_list{$cui};
	    my @vector = @{$vectors->{$vector_index}};
	    my $cos_sim = &CosSim(\@centroid, \@vector);
	    if ($cos_sim > $max_cos_sim){
		$max_cos_sim = $cos_sim; # update maximum cosine similarity value
		$centroid_vector_index = $vector_index;	# update centroid vector index
	    }
	}
	my $cui = $cui_list->{$centroid_vector_index};
	my $term = $cui_terms->{$cui};
	$cluster_names{$cluster_index} = $term;	
    }
    return \%cluster_names;
}

#Calculates cosine similarity between two vectors (aka the dot product)
# cosine similarity forumla (Vector a, Vector b):
# (a * b) / (|a|*|b|)
sub CosSim {
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


########################################
# Read Input Files
########################################

#Reads in the LBD target term fileName
# input: the LBD target term file
# output: the cui indeces, cui scores, and preferred terms
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


#Reads the vector file
# input: the list of CUIs to get vectors for and the vector fileName
# output: a list of vectors and a new list of CUIs (since not all CUIs have vectors)
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


#Gets the size of a matrix (stored as a hash of hashes)
# input: a hash of hashes
# output: the input's numRows, numCols
sub GetMatrixSize {
    my ($vectors) = (@_);
    my @matrix_size;
    $matrix_size[0] = keys %$vectors;	# number of rows
    $matrix_size[1] = @{$vectors->{0}};	# number of columns
    return \@matrix_size;
}




########################################
# VCluster Prep and Running
########################################

# converts internal format to VCluster format so it can be run
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

#runs Vcluster and outputs the clusters and the clusters tree
sub RunVCluster {
    my ($output_dir, $vcluster_dir, $v_ifile, $cl_method, $matrix_size) = (@_);

    #We want a cluster tree where each leaf node is a vector
    # thereby creating a full cluster tree. SO num_clusters = num_vectors
    my $num_vectors = @$matrix_size[0];
    my $num_clusters = $num_vectors;

    #create file names
    my $file_to_cluster = $output_dir.$v_ifile;
    my $clusters_filename = $file_to_cluster.'.clusters';
    my $tree_filename =  $file_to_cluster.'.tree';
    my $labels_filename = $file_to_cluster.'labels';

    #NOTE: the labeltree option labels with a set of features, not a single name
    #create command and perform clustering
    #my $cmd = "./$vcluster_dir/vcluster $file_to_cluster $num_clusters -clmethod=$cl_method -clustfile=$clusters_filename -showtree -cltreefile=$tree_filename";
 my $cmd = "./$vcluster_dir/vcluster $file_to_cluster $num_clusters -clmethod=$cl_method -clustfile=$clusters_filename -fulltree -treefile=$tree_filename";
    print "$cmd\n";
    my $cluster_cmd = `$cmd`;

    return $clusters_filename, $tree_filename;
}













########################################
# Get and Check Input Parameters
########################################
sub GetArgs {
    my $vcluster_dir = $ARGV[0];
    my $lbd_file = $ARGV[1];
    my $vector_file = $ARGV[2];
    my $cl_method = $ARGV[3];
    my $out_file = $ARGV[4];
    &CheckArgv($vcluster_dir, $lbd_file, $vector_file, $cl_method, $out_file);
    return $vcluster_dir, $lbd_file, $vector_file, $cl_method, $out_file;
}

sub CheckArgv {
    my ($vcluster_dir, $lbd_file, $vector_file, $cl_method, $out_file) = (@_);
    &CheckNumArgs();
    &CheckVcluster($vcluster_dir);
    &CheckCLMethod($cl_method);
    &CheckFileErr($lbd_file);
    &CheckFileErr($vector_file);
    &CheckCreateFileErr($out_file);
}

sub CheckNumArgs {
    if ($#ARGV < 4) {
	print "Incorrect number of arguments. Program requires 5 arguments to run.\n";
	print "Usage: perl DiscoveryReplication.pl [vcluster_dir] [lbd_file] [vector_file] [cl_method] [out_file]\n";
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

sub CheckCreateFileErr {
    my ($file) = (@_);
    open OUT, ">$file" or die("ERROR: cannot create file: $file\n");
    close OUT;
}

sub CreateDir {
    # checks if dir exists. if not, creates dir.
    my ($dir) = (@_);
    if (! -e $dir || ! -d $dir){
	mkdir $dir;
    }
}
