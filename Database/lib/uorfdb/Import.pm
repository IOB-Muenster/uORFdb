#!/usr/bin/env perl


# AUTHOR

# Felix Manske, felix.manske@uni-muenster.de


# COPYRIGHT

# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO,  PROCUREMENT  OF  SUBSTITUTE GOODS  OR  SERVICES;
# LOSS  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
# OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


package Import;
use Exporter;
use Encode;
use URI::Escape;
use Data::Dumper;
use List::Util qw(first);


@ISA = ('Exporter');
@EXPORT      = ();
@EXPORT_OK   = qw(readFile add2table getKey insertHash insertRowwise deleteFromTable getWhoWhen);


# Remove any type of whitespace on both ends of all strings in an array.
sub trimArray {
	my @array = @{$_[0]};
	
	foreach my $elem (@array) {
		$elem=~  s/^\s+|\s+$//g;
	}
	return \@array;
}

#
# ===============================================================================
# Basic functions to initialize the database
# ===============================================================================
#
use strict;
use warnings;
use Digest::MD5 qw(md5_hex);
use Data::Dumper;
use CGI qw(remote_addr);


# Reads in file and returns a hash with each column as the key
# and the row values for that column as an array reference.
# Columns can be (in/ex)cluded based on the index starting from 0.
# These must not be specified in the columns array.
# Speficy whether you want to include or exclude.
sub readFile {
	my $file = $_[0];
	my @columns = @{$_[1]};
		
	my $command = $_[3];
	my $keepWhitespaces = "";
	$keepWhitespaces = $_[6] if defined $_[6];
	my @uniqueIdxs = ();
	@uniqueIdxs = @{$_[4]} if defined $_[4];
	
	# If fields with these indices are empty,
	# drop the whole line
	my @notEmptyIdxs = ();
	@notEmptyIdxs = @{$_[5]} if defined $_[5];
	my $specFieldsC = 0;
	
	my $unique = "";
	
	my $insertC = 0;
	
	my @ignore = ();
	my @include = ();
	my %md5s = ();
	my $md5 = "";
	
	my $dups = 0;
	my $empty = 0;
	
	if ($command eq "ignore") {
		# Positions to ignore in each line
		@ignore = @{$_[2]};
		@ignore = sort{ $a<=>$b } @ignore;
	}
	elsif ($command eq "include") {
		# Positions to include in each line
		# By ommiting sort, user can change pos of elements
		@include = @{$_[2]};
	}
	else {
		die "ERROR: Don't know what to do with positions. Ignore/Include?\n"
	}
	# Init hash with array refs as values
	my %data = map{ $_ => [] } @columns;
	
	
	open(FILE, "<", $file) or die "Could not open $file";
	while(<FILE>) {
		
		next if ($_ =~ m/^#/);
		
		chomp($_);
		
		# Drop empty lines
		if ($_ =~ m/^\W*$/) {
			$empty++;
			next
		}
		
		my $i = 0;

		# Don't strip empty fields at the end		
		my @line = split("\t", $_, -1);
		
		# Remove lines with empty fields in a set of special fields
		if (grep {$_ =~ m/^\s*$/} @line[@notEmptyIdxs]) {
			$specFieldsC++;
			next;
		}
		
		# Unique constraint on that column
		$unique = join("#", @line[@uniqueIdxs]) if (@uniqueIdxs);
		
		# Ignore entries
		if (@ignore){
			foreach my $ignore (@ignore) {
				# Target indices must be reduced, if positions have been
				# removed before
				my $remPos = $ignore - $i;
												
				splice(@line, $remPos, 1);
				$i++;
			}
		}
		# Include only selected entries.
		elsif (@include) {
			@line = @line[@include];
		}
		
		# Remove any trailing/leading whitespace for all array elements
		@line = @{trimArray(\@line)} if (not $keepWhitespaces eq "keep_whitespaces");
				
		# Drop duplicates		
		$md5 = md5_hex(@line);
		
		if (not @uniqueIdxs) {
			if (exists $md5s{$md5}) {
				$dups++;
				next;
			}
			else {
				$md5s{$md5} = undef;
			}
		}
		else {
			# The item with unique constraint already appeard...
			if (exists $md5s{$unique}) {
				#... and the extracted values are identical --> OK, use only once
				if ($md5s{$unique}->[0] eq $md5) {
					$dups++;
					next;
				}
				#... but in another context --> error
				else {
					my $uniqPrint = $unique =~ s/#/, /gr;
					my $idx = join(", ", @uniqueIdxs);
					print "ERROR: Unique constraint on $uniqPrint (column(s) $idx): Ambigous.\n";
					print "WAS:\t$md5s{$unique}->[1]\n";
					print "IS:\t".join("\t", @line)."\n";
					die;
				}
			}
			else {
				$md5s{$unique} = [$md5, join("\t", @line)];
			}
		}
		
		# Save lines like this col1 =>[row1_col1, ... rown_col1], ...		
		for (my $i = 0; $i <= $#columns; $i++ ) {
			
			# DBI should interprete empty strings as NULL.
			my $insert = undef;
			if ($line[$i] !~ m/^$/) {
				$insert = Encode::decode('UTF-8', $line[$i]) #uri_unescape($line[$i])
			}
			push(@{$data{$columns[$i]}}, $insert);
		}
		
		$insertC ++;
	}
	
	print "WARNING: Dropped $dups duplicates\n" if ($dups != 0);
	print "WARNING: Found $empty empty lines\n" if ($empty != 0);
	print "WARNING: Removed $specFieldsC lines with empty special fields.\n" if ($specFieldsC != 0);
	return \%data;
}


# Add data from hash to table hash without overwriting.
sub add2table {
	my $tableR = $_[0];
	my %data = %{$_[1]};
	
	foreach my $key (keys(%data)) {
		if (exists $tableR->{$key}) {
			print "WARNING: Key $key already exists in table\n";
			next;
		}
		else {
			$tableR->{$key} = $data{$key}
		}
	}
	
	return $tableR;
}


# Delete rows by index from table hash.
# Useful to delete rows with null values.
sub deleteFromTable {
	
	my $tableR = $_[0];
	my @remIdxs = @{$_[1]};
	
	# Make sure that remIdxs are unique
	my %rems = map {$_ => undef} @remIdxs;
	@remIdxs = keys(%rems);
	
	@remIdxs = sort { $b <=> $a } @remIdxs;
	
	foreach my  $key (keys(%{$tableR})) {
		
		# All values should be deleted from array
		if (@remIdxs == @{$tableR->{$key}}) {
			delete $tableR->{$key};
			next;
		}
		
		# Only some values should be deleted from array
		foreach my $remIdx (@remIdxs) {
			splice (@{$tableR->{$key}}, $remIdx, 1)
		}
	}
	
	print "INFO: Deleted ".@remIdxs." records\n";
}


# Get foreign keys to connect the tables
# One query value will only return one match by default.
# An array reference of all matches can be returned by setting
# $mode to 'all'
sub getKey {
		
	my %queries = %{$_[0]};
	my $tcol = $_[1];
	my $ttable = $_[2];
	my $dbh = $_[3];
	my $mode = $_[4] || 'first-strict';
	
	# An arrayR for column names of which at least one must be defined
	my $modeOrR = $_[5];
	
	my $keepUndef = $_[6] || 0;
	
	my @keys = ();
	my @remIdxs = ();
	my %md5s = ();
	my @modeOrIdxs = ();
	
	my $unmatched = 0;
	my $missingData = 0;
	
	# Check if all queries have the same array length = row number
	my $arrayLen = undef;
	
	my @cols = keys(%queries);
	
	for (my $i = 0; $i <=$#cols; $i ++) {
		
		my $col = $cols[$i];
		
		if (defined $arrayLen) {
			if ($arrayLen != @{$queries{$col}}) {
				die "ERROR: Queries must have the same row number"
			}
		}
		else {
			$arrayLen = @{$queries{$col}};
		}
		
		# Get column indices for "mode_or"
		if (defined $modeOrR) {
			foreach my $elem (@{$modeOrR}) {
				if ($elem =~ m/^$col$/) {
					push (@modeOrIdxs, $i);
				}
			}
		}
	}
	
	# Execute row wise SELECT from column wise arrays
	my @values = values (%queries);
	
	for (my $i = 0; $i <= $#{$values[0]}; $i ++) {
		
		my $bind = "";
		
		# Setting mode_or allows one or more query columns for foreign key to be empty
		# (at least one must be filled). If more than one query column is given, 
		# columns can also point to different foreign keys.
		if (defined $modeOrR) {
			foreach my $col (@cols) {
				# Get rid of hash and anything after that (artificial column name)
				$col =~ s/#.*$//;
				$bind .= $col.'= ? or ';
			}
			$bind =~ s/or $//;
		}
		else {
			foreach my $col (@cols) {
				# Get rid of hash and anything after that (artificial column name)
				$col =~ s/#.*$//;
				$bind .= $col.'= ? and ';
			}
			$bind =~ s/and $//;
		}
		
		my $sth = $dbh->prepare("SELECT $tcol FROM $ttable WHERE $bind;");
		
				
		my @row = map {$_->[$i]}@values;
		
		
		# MD5 can only be safely used for rows without undef values
		# At the same time, an undefined $md5 indicates that the 
		# sth needs to be prepared to allow for 'is null'.
		my $md5 = undef;
		my @undefs = grep { not defined $_ }@row;
		
		if (@undefs) {
			# In strict mode, undef values are not allowed for key lookup
			# However, if column names are supplied via mode_or,
			# Only one or more of the supplied columns must be defined.
			# Regardless of the mode, at least one value must be defined.
			 if (@undefs == @row) {	
			 	$missingData++;
			 	push(@keys, undef);
				push(@remIdxs, $i);
			 	next;
			 }
			 elsif ($mode =~ m/strict/ and @modeOrIdxs) {			 	
			 	# Get undefs from custom "mode_or"
			 	my @modeUndefs = grep { not defined $_ }@row[@modeOrIdxs];
			 	
			 	# At least one of the "mode_or" entries must be defined
			 	# Still, all other entries need to be defined.
			 	if (@modeOrIdxs == @modeUndefs or @undefs != @modeUndefs) {
			 		
			 		$missingData++;
			 		push(@keys, undef);
					push(@remIdxs, $i);
			 		next;
			 	}
			 }
			 elsif ($mode =~ m/strict/) {			 	
			 	$missingData++;
			 	push(@keys, undef);
				push(@remIdxs, $i);
			 	next;
			 }
		}
		else {
			# Local encode to handle UTF-8 chars and md5
 			$md5 = md5_hex(map{Encode::encode('UTF-8',$_)}@row);
		}
		
		
		# Don't query the same data twice. Use old data.
		if (defined $md5 and exists $md5s{$md5}) {
			
			my $key = $keys[$md5s{$md5}];
			push(@keys, $key);
			
			if (not defined $key) {
				push(@remIdxs, $i);
				$unmatched++;
			}
		}
		
		else {
			
			if (not defined $md5) {
			
				my @tmp = ();
				$bind = "";
				
				for (my $j = 0; $j <=$#row; $j++) {
					
					if (defined $row[$j]) {
						
						$bind .= "$cols[$j]=? and ";
						push (@tmp, $row[$j]);
					}
					else {
						$bind .= "$cols[$j] is null and ";
						
					}
				}
				$bind =~ s/ and $//;
				@row = @tmp;
				
				$sth = $dbh->prepare("SELECT $tcol FROM $ttable WHERE $bind;");
			}
			
			$sth->execute (@row);
			
			my @selects = @{$sth->fetchall_arrayref};
			
			# Other modes than first expect array refs as values
			if ($mode !~ m/first/) {
					
				# Save only one position per row-md5
				$md5s{$md5} = $i if (defined $md5);
				
				my @res = map{$_->[0]} @selects;
				
				# No matching keys
				if (not @res) {
					
					print "DEBUG: Not found ".join("\t", @row)."\n";
					
					$unmatched++;
					push(@remIdxs, $i);
					push (@keys, undef);
				}
				else {
					push (@keys, \@res);
				}
				
				next;
			}
			
			# Number of matching rows per select
			if (@selects > 1) {
				
				print "ERROR: Ambigous query\n";
				print $bind."\n";
				print Dumper(\@row);
				die;
			}
			else {
				
				# Save only one position per row-md5
				$md5s{$md5} = $i if (defined $md5);
				
				# Mark rows with foreign key = null for removal
				if (not $selects[0]->[0]) {
					$unmatched++;
					push(@remIdxs, $i);
					push (@keys, undef);
				}
				else {
					push (@keys, $selects[0]->[0])
				}
			}
		}
	}
	# Optionally keep foreign keys with Null values.
	@remIdxs = () if ($keepUndef == 1);
	
	print "WARNING: $unmatched unmatched records\n" if ($unmatched != 0);
	print "WARNING: $missingData records with insufficient data\n" if ($missingData != 0);	
	return \@keys, \@remIdxs;
}


# Takes a hash and inserts it into a table.
# Keys must be named as column names in table.
# Values are arrays with data for that column
# over all rows
# c1 => [r1, r2,.., rN], ... cN => [r1, r2,.., rN]
sub insertHash {
	my %hash = %{$_[0]};
	my $table = $_[1];
	my $dbh = $_[2];
	my @conflCols = @{$_[3]};
	my $conflBehav = $_[4];
	
	my $bind = join(',', ('?') x keys(%hash));
	my $columns = join(',', keys(%hash));
	
	my $sth = "";
	# Skip or update conflicting entries between table hash and database
	if (@conflCols and defined $conflBehav) {
		if ($conflBehav eq "ignore") {
			$sth = $dbh->prepare("INSERT INTO $table (" . $columns . ") VALUES (" . $bind . ") ON CONFLICT (" .
			join(",", @conflCols) . ") DO NOTHING;");
		}
		elsif ($conflBehav eq "update") {
			# Prepare the update statement. Keyword EXCLUDED refers to values in VALUES()
			my @updates = map{$_ . " = EXCLUDED." . $_}keys(%hash);
			my $update = join(", ", @updates);
			
			$sth = $dbh->prepare("INSERT INTO $table (" . $columns . ") VALUES (" . $bind . ") ON CONFLICT (" .
			join(",", @conflCols) . ") DO UPDATE SET $update;");
		}
		else {
			die "ERROR: Unknown behaviour for conflicting entries ->$conflBehav<-";
		}
	}
	# Default: Normal insert
	else {
		$sth = $dbh->prepare("INSERT INTO $table (" . $columns . ") VALUES (" . $bind . ");");
	}
	
	my $tuples = $sth->execute_array({ ArrayTupleStatus => \my @tuple_status },
		values(%hash));
	
	# $tuples is undef if any tuple execution failed.
	if ($tuples) {
		my $conflict = "";
		if (@conflCols and defined $conflBehav) {
			$conflict = " Chose $conflBehav for potential conflicts in " . join("; ", @conflCols) . "."
		}
		
    	print "INFO: Successfully inserted $tuples records.$conflict\n";
	}
	else{
		print "ERROR: Failed to insert hash into $table\n";
		die;
	}
}


# Takes a hash and inserts it into a table.
# Keys must be named as column names in table.
# Values are arrays with data for that column
# over all rows. Compared to insertHash,
# insertRowwise can handle fields that contain array
# references.
# c1 => [r1, r2,.., rN], ... cN => [r1, r2,.., rN]
sub insertRowwise {
	my %hash = %{$_[0]};
	my $table = $_[1];
	my $dbh = $_[2];
	my %refCols = %{$_[3]};
	my @conflCols = @{$_[4]};
	my $conflBehav = $_[5];
	
	my $insertC = 0;
	
	my $bind = join(',', ('?') x keys(%hash));
	
	my @columns = keys(%hash);
	my $cols = join (", " ,@columns);

	# Error handling: Check that number of rows is the same for all columns
	my $rowNumb = 0;
	my $lastRowNumb = undef;
	
	foreach my $col (@columns) {
		$rowNumb = @{$hash{$col}}-1;
		die "ERROR: Column number not equal in $col" if (defined $lastRowNumb and $rowNumb != $lastRowNumb);
		$lastRowNumb = $rowNumb;
	}
	
	# Execute insert for one row at a time
	# Push the returned id for the new entry to an array
	my @res =();
	
	for (my $i = 0; $i <= $rowNumb; $i++) {
		
		my @row = map {$hash{$_}->[$i]}keys(%hash);
		
		# Some columns may contain array references which need to be prepared.
		if (keys(%refCols)) {
			
			my $i = 0;
			my $rowNumb = 0;
			my $maxRowNumb = undef;
			
			# Perform as many insert statements per row, as I have values in the array references
			while (0==0) {
				
				my @tmp = ();
				
				# Loop over the columns and get the values. Check if array refs have the same number of
				# elements (only once)
				for (my $j = 0; $j <= $#columns; $j++) {
					
					my $value = $row[$j];
					
					# It's a column with an array ref
					if (exists $refCols{$columns[$j]}) {	
						if ($i == 0) {
							$rowNumb = @{$value};
							die "ERROR: Array ref length not equal. Row: " . Dumper(\@row)
								. "\n Columns: " . $cols if (defined $maxRowNumb and $rowNumb != $maxRowNumb);
							$maxRowNumb = $rowNumb;
						}
						$value = $row[$j]->[$i];
					}
					push (@tmp, $value);
				}
				
				last if ($i > $maxRowNumb-1);
				
				my $sth = "";
				
				# Skip or update conflicting entries between table hash and database
				if (@conflCols and defined $conflBehav) {
					if ($conflBehav eq "ignore") {
						$sth = $dbh->prepare("INSERT INTO $table (" . $cols . ") VALUES (" . $bind . ") ON CONFLICT (" .
						join(",", @conflCols) . ") DO NOTHING;");
					}
					elsif ($conflBehav eq "update") {
						# Prepare the update statement. Keyword EXCLUDED refers to values in VALUES()
						my @updates = map{$_ . " = EXCLUDED." . $_}@columns;
						my $update = join(", ", @updates);
						
						$sth = $dbh->prepare("INSERT INTO $table (" . $cols . ") VALUES (" . $bind . ") ON CONFLICT (" .
						join(",", @conflCols) . ") DO UPDATE SET $update;");
					}
					else {
						die "ERROR: Unknown behaviour for conflicting entries ->$conflBehav<-";
					}
				}
				# Default: Normal insert
				else {
					$sth = $dbh->prepare("INSERT INTO $table (".$cols.") VALUES (".$bind.");");
				}
				
				$sth->execute(@tmp);
				
				$i++;
				$insertC++;
			}			
		}
		else {
			
			my $sth = "";
			# Skip or update conflicting entries between table hash and database
			if (@conflCols and defined $conflBehav) {
				if ($conflBehav eq "ignore") {
					$sth = $dbh->prepare("INSERT INTO $table (" . $cols . ") VALUES (" . $bind . ") ON CONFLICT (" .
					join(",", @conflCols) . ") DO NOTHING;");
				}
				elsif ($conflBehav eq "update") {
					# Prepare the update statement. Keyword EXCLUDED refers to values in VALUES()
					my @updates = map{$_ . " = EXCLUDED." . $_}@columns;
					my $update = join(", ", @updates);
					
					$sth = $dbh->prepare("INSERT INTO $table (" . $cols . ") VALUES (" . $bind . ") ON CONFLICT (" .
					join(",", @conflCols) . ") DO UPDATE SET $update;");
				}
				else {
					die "ERROR: Unknown behaviour for conflicting entries ->$conflBehav<-";
				}
			}
			# Default: Normal insert
			else {
				$sth = $dbh->prepare("INSERT INTO $table (".$cols.") VALUES (".$bind.");");
			}
				
			$sth->execute(@row);
			
			$insertC++;
		}
	}
	my $conflict = "";
	if (@conflCols and defined $conflBehav) {
		$conflict = " Chose $conflBehav for potential conflicts in " . join("; ", @conflCols) . "."
	}
	
	print "INFO: Successfully inserted $insertC records.$conflict\n";
}

1;