#!/usr/local/bin/perl -w

use DBI;
use Getopt::Long;
use FindBin qw( $Bin );

use lib "$Bin/../../perl";
use SBEAMS::Connection::Settings qw( $DBCONFIG $DBINSTANCE );
use strict;

use constant COMMIT => 20;

$|++; # don't buffer output

{ # MAIN

  # If no args are given, print usage and exit
  if (! @ARGV) {
    printUsage();
    exit 0;
  }

  my $args = processArgs();
  my $cmds = parseFile ( $args );
  if ( $args->{test_mode} ) {
    printCommands( $args, $cmds );
    exit;
  }
  my $dbh = dbConnect( $args );

 if ( $args->{query_mode} ) {
    printResults( $args, $dbh, $cmds );
  } else {
    insertRecords( $args, $dbh, $cmds );
  }
  $dbh->disconnect();
  exit 0;
}


sub parseFile {
  my $args = shift;
  my @cmds;
  my $cmd = '';

  # Manual queries take precedence...
  if ( $args->{manual_mode} ) {
    # This is gonna be querymode.
    $args->{query} ||= 1;
    print "Running manual query\n";
    return( [ $args->{manual} ] );
  }

  # Otherwise, parse the file as usual.
  print "Parsing command file $args->{sfile}\n" if $args->{verbose};
  open( FIL, $args->{sfile} ) || die "Unable to open file $args->{sfile}";
  while ( my $line = <FIL> ) {
    #chomp $line;
    next if $line =~ /^\s*$/;

    #### If the user opts to ignore Audit Trail FOREIGN KEYS, stop when found
    if ($args->{no_audit_constraints} && $line =~ /Audit trail FOREIGN KEYS/) {
      last;
    }

    if ( $args->{delimiter} eq 'GO' ) {
      if ( $line =~ /^GO\s*$/i ) {
        push @cmds, $cmd;
        $cmd = '';
      } else {
        $cmd .= "$line\n";
      }
    } else {
      if ( $line =~ /;\s*$/ ) {
        $cmd .= $line;
        push @cmds, $cmd;
        $cmd = '';
      } else {
        $cmd .= $line;
      }
    }
  }
  push @cmds, $cmd if $cmd;  # Leftovers
  return( \@cmds );
}

sub printCommands {
  my $args = shift;
  my $sql = shift;
  my $cnt = 0;
  foreach my $sql ( @$sql ) {
    $cnt++;
    print "--Stmt $cnt\n";
    print "$sql\n";
  }
}

sub printResults {
  my $args = shift;
  my $dbh = shift;
  my $sql = shift;
  foreach my $sql ( @$sql ) {
    if ( $sql !~ / UPDATE | DROP | INSERT | DELETE | TRUNCATE /gi ) {
      print "$sql\n";  # Query mode implies a certain amount of verbosity!
      my $sth = $dbh->prepare( $sql );
      $sth->execute();
      my $firstrow = 1;
      while ( my $row = $sth->fetchrow_hashref() ) {
        if ( $firstrow ) {
          print join( "\t", keys( %{$row} ) . "\n" );
          $firstrow = 0;
        }
        print join( "\t", values( %{$row} ) . "\n" );
      }
    } else {
      print "Query mode is read-only, not running $sql";
    }
  }
}

sub insertRecords {
  my $args = shift;
  my $dbh = shift;
  my $cmds = shift;
  
  my $cnt;
  foreach my $cmd ( @$cmds ) 
  {
      $cnt++;
      print "$cmd\n" if $args->{verbose};

      if ($args->{verbose})
      {
          eval 
          {
              $dbh->do( $cmd );
              $dbh->commit();
          };
          if( $@ ) 
          {
              $dbh->rollback();
          }
      } else 
      {
          $dbh->do( $cmd );
      }
  }
}

sub dbConnect {
  my $args = shift;

  # Get db driver from configuration file
  my $DB_SERVER = $DBCONFIG->{$DBINSTANCE}->{DB_SERVER};
  my $DB_DATABASE = $args->{database} || $DBCONFIG->{$DBINSTANCE}->{DB_DATABASE};
  my $cstring = eval "\"$DBCONFIG->{$DBINSTANCE}->{DB_DRIVER}\"";

  my $dbh = DBI->connect( $cstring, $args->{user}, $args->{pass} ) || die ('couldn\'t connect' );
#  $dbh->{AutoCommit} = 0;
  $dbh->{RaiseError} = ( $args->{ignore_errors} ) ? 0 : 1;

  if ($args->{fudge})
  {
      $dbh->{RaiseError} = 0; ## no errors raised
      $dbh->{AutoCommit} = 0; ## no auto commit
      $dbh->{PrintError} = 0; ## no errors printed
  }

  print "Connected to database successfully\n" if $args->{verbose};
  return $dbh;
}

sub processArgs {
  my %args;
  unless( GetOptions ( \%args, 'pass=s', 'user=s', 'verbose', 'sfile=s',
                      'delimiter=s', 'ignore_errors', 'manual:s',
                      'no_audit_constraints', 'database=s', 'query_mode', 
                      'test_mode', 'fudge' ) ) {
  printUsage("Error with options, please check usage:");
  }

  for ( qw( user ) ) {
    if ( !defined $args{$_} ) {
    printUsage( "Missing required parameter: $_" ) unless $args{test_mode}; 
    }
  }
  unless( $args{manual} || $args{sfile} ) {
    printUsage( "Must specify either a sql file or a manual query" );
  }

  # User declined to enter a password, prompt for one
  while ( !$args{pass} ) {
    if ($args{test_mode}) {
      $args{pass} = 'testing';
      next;
    }
    print "Enter password, followed by [Enter] (cntl-C to quit):\t";
    $|++;
    system("stty -echo");
    my $pass = <>;
    system("stty echo");
    chomp $pass;
    print "\n";
    if ( $pass ) {
      #( my $err = $pass ) =~ s/[\w\#]//g;
      #die "Illegal characters in password: $err \n" if $err;
      $args{pass} = $pass;
    } else {
      print "No input received.\n";
      $|++;
    }
  }
  
  # Delimiter will either be semicolon or GO
  $args{delimiter} = ( !$args{delimiter} ) ? ';' :
                     ( uc($args{delimiter}) eq 'GO' ) ? 'GO' : ';'; 

  return \%args;
}

sub printUsage {
  my $err = shift || '';
  print( <<"  EOU" );
   $err
   
   Usage: $0 -u username -s sfile [ -p password ]

   -u --user xxxx     username to authenticate to the db
   -p --pass xxxx     password to authenticate to the db.  Will be prompted
                      if value is ommitted.
   -s --sfile xxxx    SQL file which defines table and columns etc
   -v --verbose       verbose output
   -i --ignore_errs   Ignore SQL errors and continue
   -d --delimiter xx  Delimter for splitting file, semicolon (default) or GO.
   -q --query_mode    Run (SELECT) query(s) and print results.
   -t --test_mode     Parse file and simply print out each statement that would have been executed. 
   -m --manual        SELECT query provided explicitly, obviates the need for
                      a SQL file.
   -n --no_audit_constraints  If set, then the Audit Trail FOREIGN KEYS are skipped
      --database      Specify a database to initially connect to besides the default
   -f --fudge         Current design of tables is such that new installations may
                      have modules with tables in common.  Create table and constraint
                      statements no that case leads to a flood of errors, which
                      are inherent to the design.  Silencing those errors here is
                      a fudge until table design changes.

  EOU
  exit;
}
