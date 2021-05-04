package SBEAMS::PeptideAtlas::HTMLPrinter;

###############################################################################
# Program     : SBEAMS::PeptideAtlas::HTMLPrinter
# Author      : Eric Deutsch <edeutsch@systemsbiology.org>
# $Id$
#
# Description : This is part of the SBEAMS::WebInterface module which handles
#               standardized parts of generating HTML.
#
#		This really begs to get a lot more object oriented such that
#		there are several different contexts under which the a user
#		can be in, and the header, button bar, etc. vary by context
###############################################################################
use strict;
use LWP::UserAgent;
use HTTP::Request;
use Data::Dumper;
use LWP::Simple;

use vars qw($sbeams $current_contact_id $current_username $q
             $current_work_group_id $current_work_group_name
             $current_project_id $current_project_name $current_user_context_id);
use CGI::Carp qw(fatalsToBrowser croak);
use SBEAMS::Connection qw($log $q);
use SBEAMS::Connection::DBConnector;
use SBEAMS::Connection::Settings;
use SBEAMS::Connection::TableInfo;
use SBEAMS::Connection::DataTable;
use SBEAMS::Connection::GoogleVisualization;

use SBEAMS::PeptideAtlas::Settings;
use SBEAMS::PeptideAtlas::TableInfo;
use SBEAMS::PeptideAtlas::Tables;

use SBEAMS::Proteomics::Tables;

###############################################################################
# printPageHeader
###############################################################################
sub printPageHeader {
  my $self = shift;
  $self->display_page_header(@_);
}


###############################################################################
# display_page_header
###############################################################################
sub display_page_header {
  my $self = shift;

  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  my $sbeams = $self->getSBEAMS();
  my $project_id = $args{'project_id'} || $sbeams->getCurrent_project_id();
  my $output_mode = $sbeams->output_mode();
  my $http_header = $sbeams->get_http_header();

  #### If the output mode is interactive text, display text header
  if ($output_mode eq 'interactive') {
    $sbeams->printTextHeader();
    return;
  }
  #### If the output mode is interactive text, display text header
  if ($output_mode =~ 'xml') {
    my $xml_header = $sbeams->get_http_header( mode => 'xml', filename => 'peptide_export.xml' );
    print STDERR $xml_header;
    print $xml_header;
    return;
  }#elsif ($output_mode =~ 'tsv'){
  #  my $text_header = $sbeams->get_http_header( mode => 'tsv', filename => 'table.tsv' );
  #  print STDERR $text_header;
  #  print $text_header;
  #  return 
  #}

  #### If the output mode is not html, then we may not want a header here
  if ($output_mode ne 'html') {

    # Caller may want header printing to be handled here
    if ( $args{force_header} ) { 

      # Print http header
      print $http_header if $sbeams->invocation_mode() eq 'http';
    }

    return;
  }

  #### Obtain main SBEAMS object and use its http_header
  $sbeams = $self->getSBEAMS();

  if( $sbeams->isGuestUser() || $sbeams->getCurrent_username() =~ /^reviewer/i ) {
      $self->displayGuestPageHeader( @_ );
      return;
  } elsif ( $CONFIG_SETTING{PA_USER_SKIN} ) {
    $current_username = $sbeams->getCurrent_username();
    for my $skin ( split( ",", $CONFIG_SETTING{PA_USER_SKIN} ) ) {
#      $log->debug( "skin is $skin" );
      my ( $name, $value ) = split( /::::/, $skin, -1 );
#      $log->debug( "name is $name, val is $value" );
      if ( $name eq $current_username ) {
#        $log->debug( "CUSTOM! $CONFIG_SETTING{PA_USER_SKIN}" );
        $self->displayGuestPageHeader( @_, uri => $value );
        return;
      } else {
        $log->debug( "$name does not equal $current_username" );
      }
    }
  }
  $self->displayStandardPageHeader(@_);
}


###############################################################################
# displayInternalResearcherPageHeader
###############################################################################
sub displayInternalResearcherPageHeader {
  my $self = shift;
  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  my $LOGOUT_URI = "$CGI_BASE_DIR/logout.cgi";

  my $LOGOUT_LINK = qq~<A HREF="$LOGOUT_URI" class="Nav_link">LOGOUT</A>~;


  #### Obtain main SBEAMS object and use its http_header
  my $sbeams = $self->getSBEAMS();
  my $http_header = $sbeams->get_http_header();
  use LWP::UserAgent;
  use HTTP::Request;
  my $ua = LWP::UserAgent->new();
  my $skinLink = 'http://www.peptideatlas.org';
  my $response = $ua->request( HTTP::Request->new( GET => "$skinLink/.index.dbbrowse.php" ) );
  my @page = split( "\n", $response->content() );
  my $skin = '';
  my $cnt=0;
  for ( @page ) {

    $cnt++;
    $_ =~ s/\<\!-- LOGIN_LINK --\>/$LOGOUT_LINK/;
    last if $_ =~ /--- Main Page Content ---/;
    $skin .= $_;
  }

  $skin =~ s/images\//sbeams\/images\//gm;
#  $skin =~ s/\/images\//\/sbeams\/images\//gm;
  $self->{'_external_footer'}=join("\n", @page[$cnt..$#page]);
  $self->{'_external_footer'} =~ s/images\//sbeams\/images\//gm;
 
  print "$http_header\n\n";
  print <<"  END_PAGE";
  $skin
  END_PAGE

  $self->printJavascriptFunctions();
  }


###############################################################################
# displayGuestPageHeader
###############################################################################
sub displayGuestPageHeader {
 	my $self = shift;
  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  my $LOGIN_URI = "$SERVER_BASE_DIR$ENV{REQUEST_URI}";

  if ($LOGIN_URI =~ /\?/) {
    $LOGIN_URI .= ";force_login=yes";
  } else {
    $LOGIN_URI .= "?force_login=yes";
  }
  my $LOGIN_LINK = qq~<A HREF="$LOGIN_URI" class="Nav_link">LOGIN</A>~;

  my $doctype = '<!DOCTYPE HTML>' . "\n";
  $doctype = '' unless $args{show_doctype}; 


	my $sbeams = $self->getSBEAMS();
  my $message = $sbeams->get_page_message( q => $q );
 	$current_username ||= $sbeams->getCurrent_username();


	my $cswitcher = '';
	if ( -e "$PHYSICAL_BASE_DIR/lib/perl/SBEAMS/PeptideAtlas/ContextWidget.pm" ) {
		require SBEAMS::PeptideAtlas::ContextWidget;
		my $cwidget = SBEAMS::PeptideAtlas::ContextWidget->new();
		$cswitcher = $cwidget->getContextSwitcher( username => $current_username,
		                                           cookie_path => $HTML_BASE_DIR );
	} else {
		$log->debug(  "$PHYSICAL_BASE_DIR/lib/perl/SBEAMS/PeptideAtlas/ContextWidget.pm doesn't exist" ) 
	}


  # Use http_header from main SBEAMS object
  my $http_header = $sbeams->get_http_header();

  my $js =  "<SCRIPT LANGUAGE=javascript SRC=\"/sbeams/usr/javascript/sorttable.js\"></SCRIPT>";


  my $ua = LWP::UserAgent->new();
  my $skinLink = $args{uri} || 'http://www.peptideatlas.org/.index.dbbrowse.php';
  my $resource = $sbeams->getSessionAttribute( key => 'PA_resource' ) || '';
  if ( $resource eq 'SRMAtlas' ) {
    $skinLink = 'http://www.srmatlas.org/.index.dbbrowse.php';
  } elsif ( $resource eq 'DIAAtlas' ) {
#    $skinLink = 'http://www.srmatlas.org/.index.dbbrowse.php';
    $skinLink = 'http://www.swathatlas.org/.index.dbbrowse.php';
  }
  my $response = $ua->request( HTTP::Request->new( GET => "$skinLink" ) );
  my @page = split( "\n", $response->content() );
	if ( $args{show_doctype} ) {
    my $first = shift @page;
		if ( $first !~ /doctype/i ) {
		  unshift @page, $first;
		}
		unshift @page, $doctype;
	}

  my $skin = '';
  my $cnt=0;
  my $init = ( $args{init_tooltip} ) ? $self->init_pa_tooltip() : '';
  my $css_info = $sbeams->printStyleSheet( module_only => 1 );
  my $loadscript = "$args{onload};" if $args{onload};
  $loadscript .= 'sortables_init();' if $args{sort_tables};

  $LOGIN_LINK .= "<BR><BR><BR>\n$cswitcher<BR>\n";
  for ( @page ) {
    $cnt++;

    # Login link originates in Peptide/SRM/SWATH Atlas. This section triages link
    # to work for current dev instance
    if ( $_ =~ /force_login=yes/ ) {
      my $url = $q->self_url();
      my $sep = ( $url =~ /\?/ ) ? ';' : '?';
      $url = $url . $sep . 'force_login=yes' unless $url =~ /force_login/;
      my $site_url = $_;
      $site_url =~ s/^(.*HREF=")[^"]+(".*$)/$1$url$2/g;
      $skin .= $site_url;
      next;
    }
#   Turned of this mechanism - was working only for Peptide Atlas, and
#   this was printing a second login link.
#    $_ =~ s/\<\!-- LOGIN_LINK --\>/$LOGIN_LINK/;

    $_ =~ s/(\<[^>]*body[^>]*\>)/$1$init$css_info/;
    if($loadscript){
      $_ =~ s/<body/$js\n<body OnLoad="$loadscript self.focus();"/;
    }
    $_ =~ s/width="680"/width="100%"/;  # resultsets are often wide...
    $_ =~ s/width="550"//;              # and IE has trouble rendering these tables
    last if $_ =~ /--- Main Page Content ---/;
    $skin .= "$_\n";
  }
 
  $self->{'_external_footer'} = join("\n", '<!--SBEAMS_PAGE_OK-->', @page[$cnt..$#page]);
  $skin =~ s#/images/#/sbeams/images/#gm;
  $self->{'_external_footer'} =~ s#/images/#/sbeams/images/#gm;
  #$skin =~ s#/images/#/dev2/sbeams/images/#gm;
  
  if ( $args{tracker_type} ) {
    if ( $args{tracker_type} eq 'srm' ) {
#      $skin =~ s/<!-- Google Analytics.*<!-- End Google Analytics -->//gms;
#      die $skin;
    } elsif ( $args{tracker_type} eq 'swath' ) {
      if ( $skinLink !~ /swath/ ) {
        $skin =~ s/<!-- Google Analytics.*<!-- End Google Analytics -->//gms;
        $self->{'_external_footer'} =~ s/<!-- Google Analytics.*<!-- End Google Analytics -->//gms;
        my $track_link = $skinLink;
        $track_link =~ s/.index.*$/includes\/tracker.inc.php/;
        my $response = $ua->request( HTTP::Request->new( GET => "$track_link" ) );
        my $tracker = $response->content();
        $skin =~ s/<head>/<head>\n$tracker/gi;
      }
    }
  }

  print "$http_header\n\n";
  print "$skin";
  print "$args{header_info}\n" if $args{header_info};

  $self->printJavascriptFunctions();
  print $message;
  }



###############################################################################
# displayStandardPageHeader
###############################################################################
sub displayStandardPageHeader {
  my $self = shift;
  my %args = @_;

  my $navigation_bar = $args{'navigation_bar'} || "YES";

  #### Obtain main SBEAMS object and use its http_header
  my $sbeams = $self->getSBEAMS();
  my $http_header = $sbeams->get_http_header();

  my $message = $sbeams->get_page_message( q => $q );

  my $doctype = '<!DOCTYPE HTML>' . "\n";
  $doctype = '' unless $args{show_doctype}; 

  print qq~$http_header
	$doctype<HTML><HEAD>
	<TITLE>$DBTITLE - $SBEAMS_PART</TITLE>
  ~;

  $self->printJavascriptFunctions();
  $self->printStyleSheet();

  my $loadscript = "$args{onload};" || '';
  $loadscript .= 'sortables_init();' if $args{sort_tables};
  my $js =  "<SCRIPT LANGUAGE=javascript SRC=\"/sbeams/usr/javascript/sorttable.js\"></SCRIPT>";
  print "$args{header_info}\n" if $args{header_info};

  #### Determine the Title bar background decoration
  my $header_bkg = "bgcolor=\"$BGCOLOR\"";
  $header_bkg = "background=\"$HTML_BASE_DIR//images/plaintop.jpg\"" if ($DBVERSION =~ /Primary/);

  print qq~
	<META http-equiv="X-UA-Compatible" content="chrome=IE8,IE=edge">
	<!--META HTTP-EQUIV="Expires" CONTENT="Fri, Jun 12 1981 08:20:00 GMT"-->
	<!--META HTTP-EQUIV="Pragma" CONTENT="no-cache"-->
	<!--META HTTP-EQUIV="Cache-Control" CONTENT="no-cache"-->
	</HEAD>

	<!-- Background white, links blue (unvisited), navy (visited), red (active) -->
  ~;
  if($loadscript){
     print qq~
       $js
	     <BODY BGCOLOR="#FFFFFF" TEXT="#000000" LINK="#0000FF" VLINK="#000080" ALINK="#FF0000" TOPMARGIN=0 LEFTMARGIN=0 OnLoad="$loadscript self.focus();">
     ~;
  }
  else{
      print qq~
              <BODY BGCOLOR="#FFFFFF" TEXT="#000000" LINK="#0000FF" VLINK="#000080" ALINK="#FF0000" TOPMARGIN=0 LEFTMARGIN=0 OnLoad="self.focus();">
     ~;
  }
  print $self->init_pa_tooltip() if $args{init_tooltip};

  print qq~
	<table border=0 width="100%" cellspacing=0 cellpadding=1>

	<!------- Header ------------------------------------------------>
	<a name="TOP"></a>
	<tr>
	  <td bgcolor="$BARCOLOR"><a href="http://db.systemsbiology.net/"><img height=64 width=64 border=0 alt="ISB DB" src="$HTML_BASE_DIR/images/dbsmlclear.gif"></a><a href="https://db.systemsbiology.net/sbeams/cgi/main.cgi"><img height=64 width=64 border=0 alt="SBEAMS" src="$HTML_BASE_DIR/images/sbeamssmlclear.gif"></a></td>
	  <td align="left" $header_bkg><H1>$DBTITLE - $SBEAMS_PART<BR>$DBVERSION</H1></td>
	</tr>

  ~;

  #print ">>>http_header=$http_header<BR>\n";

  if ($navigation_bar eq "YES") {

    print qq~
	<!------- Button Bar -------------------------------------------->
	<tr><td bgcolor="$BARCOLOR" align="left" valign="top">
	<table border=0 width="120" cellpadding=2 cellspacing=0>

	<tr><td><a href="$CGI_BASE_DIR/main.cgi">$DBTITLE Home</a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_PART/main.cgi">$SBEAMS_PART Home</a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/logout.cgi">Logout</a></td></tr>
	<tr><td>&nbsp;</td></tr>
	<tr><td>Browse Data:</td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetPeptides"><nobr>&nbsp;&nbsp;&nbsp;Browse Peptides</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/BrowseBioSequence.cgi"><nobr>&nbsp;&nbsp;&nbsp;Browse Bioseqs</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetPeptide"><nobr>&nbsp;&nbsp;&nbsp;Get Peptide Summary</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetProtein"><nobr>&nbsp;&nbsp;&nbsp;Get Protein Summary</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/showPathways"><nobr>&nbsp;&nbsp;&nbsp;Pathway Search</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/GetTransitions"><nobr>&nbsp;&nbsp;&nbsp;SRM Transitions</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/Glycopeptide/Glyco_prediction.cgi"><nobr>&nbsp;&nbsp;&nbsp;Search Glyco-Peptides</nobr></a></td></tr>
	<tr><td>&nbsp;</td></tr>
	<tr><td>Manage Tables:</td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_biosequence_set"><nobr>&nbsp;&nbsp;&nbsp;Biosequence Sets</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_atlas_build"><nobr>&nbsp;&nbsp;&nbsp;Atlas Builds</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_default_atlas_build"><nobr>&nbsp;&nbsp;&nbsp;Default Builds</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_atlas_search_batch"><nobr>&nbsp;&nbsp;&nbsp;Atlas Search Batches</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_sample"><nobr>&nbsp;&nbsp;&nbsp;Samples</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=BL_dbxref"><nobr>&nbsp;&nbsp;&nbsp;DB Xrefs</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=organism"><nobr>&nbsp;&nbsp;&nbsp;Organisms</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_publication"><nobr>&nbsp;&nbsp;&nbsp;Publications</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_elution_time_type"><nobr>&nbsp;&nbsp;&nbsp;Elution Time Type</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_instrument_type"><nobr>&nbsp;&nbsp;&nbsp;Instrument Type</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_fragmentation_type"><nobr>&nbsp;&nbsp;&nbsp;Fragmentation Type</nobr></a></td></tr>

	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_PASS_submitter"><nobr>&nbsp;&nbsp;&nbsp;PASS</nobr></a></td></tr>
  <tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_public_data_repository"><nobr>&nbsp;&nbsp;&nbsp;Public data repository</nobr></a></td></tr>
	<tr><td>&nbsp;</td></tr>
	<tr><td>Annotations:</td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_spectrum_annotation"><nobr>&nbsp;&nbsp;&nbsp;Spectra</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_peptide_annotation"><nobr>&nbsp;&nbsp;&nbsp;Peptides</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_modified_peptide_annotation"><nobr>&nbsp;&nbsp;&nbsp;Modified Peptides</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_spectrum_annotation_level"><nobr>&nbsp;&nbsp;&nbsp;Spectra Levels</nobr></a></td></tr>
	<tr><td><a href="$CGI_BASE_DIR/$SBEAMS_SUBDIR/ManageTable.cgi?TABLE_NAME=AT_transition_suitability_level"><nobr>&nbsp;&nbsp;&nbsp;Transition Levels</nobr></a></td></tr>
	</table>
	</td>

	<!-------- Main Page ------------------------------------------->
	<td valign=top>
	<table border=0 bgcolor="#ffffff" cellpadding=4>
	<tr><td>$message

    ~;
  } else {
      print qq~
	  </TABLE>
      ~;
  }

}

# 	<table border=0 width="680" bgcolor="#ffffff" cellpadding=4>


###############################################################################
# printStyleSheet
#
# Print the standard style sheet for pages.  Use a font size of 10pt if
# remote client is on Windows, else use 12pt.  This ends up making fonts
# appear the same size on Windows+IE and Linux+Netscape.  Other tweaks for
# different browsers might be appropriate.
###############################################################################
sub printStyleSheet {
    my $self = shift;

    #### Obtain main SBEAMS object and use its style sheet
    $sbeams = $self->getSBEAMS();
    $sbeams->printStyleSheet();
    print <<"    END";
    <STYLE>
     TD { background-repeat: no-repeat; }
    </STYLE>
    END

}


###############################################################################
# printJavascriptFunctions
#
# Print the standard Javascript functions that should appear at the top of
# most pages.  There probably should be some customization allowance here.
# Not sure how to design that yet.
###############################################################################
sub printJavascriptFunctions {
    my $self = shift;
    my $javascript_includes = shift;


    print qq~
	<SCRIPT LANGUAGE="JavaScript">
	<!--

	function refreshDocument() {
            //confirm( "apply_action ="+document.MainForm.apply_action.options[0].selected+"=");
            document.MainForm.apply_action_hidden.value = "REFRESH";
            document.MainForm.action.value = "REFRESH";
	    document.MainForm.submit();
	} // end refreshDocument


	function showPassed(input_field) {
            //confirm( "input_field ="+input_field+"=");
            confirm( "selected option ="+document.forms[0].slide_id.options[document.forms[0].slide_id.selectedIndex].text+"=");
	    return;
	} // end showPassed



        // -->
        </SCRIPT>
    ~;

}


###############################################################################
# printPageFooter
###############################################################################
sub printPageFooter {
  my $self = shift;
  $self->display_page_footer(@_);
}


###############################################################################
# display_page_footer
###############################################################################
sub display_page_footer {
  my $self = shift;
  my %args = @_;


  #### If the output mode is interactive text, display text header
  my $sbeams = $self->getSBEAMS();
  if ($sbeams->output_mode() eq 'interactive') {
    $sbeams->printTextHeader(%args);
    return;
  }


  #### If the output mode is not html, then we don't want a header here
  if ($sbeams->output_mode() ne 'html') {
    return;
  }


  #### Process the arguments list
  my $close_tables = $args{'close_tables'} || 'YES';
  my $display_footer = $args{'display_footer'} || 'YES';
  my $separator_bar = $args{'separator_bar'} || 'NO';

   if($self->{'_external_footer'}) {


	  # Have to fish for error file, since we are using Carp (can't pass sbeams
		# object ).
	  my $errfile = 'Error-' . getppid();
 	  my $is_error = 0;
 	  if ( $sbeams->doesSBEAMSTempFileExist( filename => $errfile ) ) {
			$sbeams->deleteSBEAMSTempFile( filename => $errfile );
			$is_error++;
		}

   	$self->{'_external_footer'} =~ s/SBEAMS_PAGE_OK/SBEAMS_PAGE_ERROR/gmi if $is_error;

   	print "$self->{'_external_footer'}\n";
     	return;
  }


  #### If closing the content tables is desired
  if ($close_tables eq 'YES') {
    print qq~
	</TD></TR></TABLE>
	</TD></TR></TABLE>
    ~;
  }


  #### If displaying a fat bar separtor is desired
  if ($separator_bar eq 'YES') {
    print "<BR><HR SIZE=5 NOSHADE><BR>\n";
  }


  #### If finishing up the page completely is desired
  if ($display_footer eq 'YES') {
    my $close_main_tables = 'NO';
    $close_main_tables = 'YES' unless ($close_tables eq 'YES');

    #### Default to the Core footer
    $sbeams->display_page_footer(display_footer=>'YES',
      close_tables=>$close_main_tables);
  }

}

#+
# @narg header_text - Text to render as HTML
# @narg header_element - XML friendly header element
# @narg anchor - text to create <A> anchor in HTML mode
# @narg bold - header bold text
# @narg list_items - ref to array of anon hashes of form key=>value for two column table
# @narg width - fixed pixel width for table
# @narg header_width - fixed pixel width for heading
# @narg key_width - fixed percentage of width for key
sub encodeFullSectionList {
  my $self = shift || die ("self not passed");
  my %args = ( width => 600,
               header_width => 900,
               bold => 1,
               key_width => 20,
               @_
             );
  
  # Default to BOLD
  unless ( $args{header_text} && $args{list_items} ) {
    $log->error( "Required parameters not supplied" );
    return '';
  }
  $args{header_element} ||= $args{header_text};

  my $buffer;
  my $sbeams = $self->getSBEAMS();
  if ( $sbeams->output_mode() =~ /html/i ) {
    $buffer = "<table width=$args{width}>\n";
    $buffer .= $self->encodeSectionHeader( %args, text => $args{header_text} );
    for my $item ( @{$args{list_items}} ) {
      $buffer .= $self->encodeSectionItem( key => $item->{key}, value => $item->{value}, key_width => $args{key_width} . "%" ) . "\n";
    } 
    $buffer .= "</table>\n";

  } elsif ( $sbeams->output_mode() =~ /xml/i ) {
    $buffer = "<$args{header_element}>\n";
    for my $item ( @{$args{list_items}} ) {
      my $key = $item->{key};
      $key =~ s/\s//g;
      $buffer .= "<li $key='$item->{value}'/>\n";
    } 
    $buffer .= "</$args{header_element}>\n";
  } else {
    $buffer = 'Sorry, yet to be implemented!';
  }
  return $buffer;
}


###############################################################################
# encodeSectionHeader
###############################################################################
sub encodeSectionHeader {
  my $METHOD = 'encodeSectionHeader';
  my $self = shift || die ("self not passed");
  my %args = @_;

  my $colspan = $args{colspan} || 2;
  
  # Default to BOLD
  $args{bold} = 1 if !defined $args{bold};

  my $text = $args{text} || '';
  $text = "<B>$text</B>" if $args{bold};

  my $link = $args{link} || '';

  $args{anchor} ||= $args{text};
  my $anchor = ( $args{anchor} ) ? "<A NAME='$args{anchor}'></A>" : '';

  my $mouseover = '';
  if ( $args{mouseover} ) {
    my $qmark = "<img src='$HTML_BASE_DIR/images/greyqmark.gif' />";
    $text = "<div title='$args{mouseover}'>$text&nbsp;$qmark</div>";
  }

#        <TR><TD colspan="2" background="$HTML_BASE_DIR/images/fade_orange_header_2.png" width="600">$link<font color="white">$anchor$text</font></TD></TR>
#        <TR><TD colspan="2" class=fade_header width="600">$link<font color="white">$anchor$text</font></TD></TR>

  my $buffer = qq~
        <TR><TD colspan="$colspan" style="background-repeat: no-repeat; background-image: url('$HTML_BASE_DIR/images/fade_orange_header_2.png')" width="600">$link<font color="white">$anchor$text</font></TD></TR>
~;

  if ( $args{LMTABS} ) {
    my $sbeams = $self->getSBEAMS();
    $args{no_toggle} ||= 0;
    $args{divname} ||= $sbeams->getRandomString( num_chars => 20 );
    $buffer = $sbeams->make_toggle_section( neutraltext => $text,
					    sticky => 1,
					    name => $args{divname},
					    tooltip => 'Show/Hide Section',
					    barlink => 1,
					    no_toggle => $args{no_toggle},
					    visible => 1,
					    content => $text
	);
  }

###        <TR><TD colspan="$colspan" $link style="background:#f3f1e4;color:#555;border-top:1px solid #b00;border-left:15px solid #b00;padding:0.5em">$anchor$text</TD></TR> e2dcc2


  return $buffer;
  return "<table>$buffer</table>";
}


###############################################################################
# encodeSectionItem
###############################################################################
sub encodeSectionItem {
  my $METHOD = 'encodeSectionItem';
  my $self = shift || die ("self not passed");
  my %args = @_;

  my $key = $args{key} || '';

  my $value = $args{value};
  $value = '' if !defined $value;

  my $desc_nowrap = ( $args{nowrap_description} ) ? 'NOWRAP' : '';

  my $url = $args{url} || '';
  my $kwid = ( $args{key_width} ) ? "WIDTH='$args{key_width}'" : '';
  my $vwid = ( $args{val_width} ) ? "WIDTH='$args{val_width}'" : '';
  
  my $tr = $args{tr_info} || ''; 

  $url =~ s/ /+/g;
  my $astart = '<b>';
  my $aend = '</b>';
  if ($url) {
    $astart = qq~<A HREF="$url" target="_blank">~;
    $aend = qq~</A>~;
  }

  my $buffer = '';
  if ($tr =~ /hoverable/) {  # new L&F?
# e2dcc2
    $buffer = qq~
        <TR $tr><TD class="key" NOWRAP $kwid>$key</TD><TD class="value" $desc_nowrap $vwid>$astart$value$aend</TD></TR>
~;
  }
  else {
    $buffer = qq~
        <TR $tr><TD NOWRAP bgcolor="cccccc" $kwid>$key</TD><TD $desc_nowrap $vwid>$astart$value$aend</TD></TR>
~;
  }


  return $buffer;

}


###############################################################################
# encodeSectionTable
###############################################################################
#
#
# @nparam  rows             Reference to array of row arrayrefs.  Required.
# @nparam  header           rows includes header info, default 0
# @nparam  width            Width of table
# @nparam  nocenter         Overrides centering of header info
# @nparam  align            Ref to array of positional info for columns -
#                           right, left, and center
# @nparam  rows_to_show     Number of rows to show initially, others are loaded
#                           into page but hidden by CSS and showable view JS.
# @nparam  nowrap           ref to array of col indexes on which to set nowrap
# @nparam  maxrows          Maximum number of rows the table can contain, 
#                           irrespective of the number displayed.
# @nparam  chg_bkg_idx      Index of field on which to trigger alternating
#                           colors.
# @nparam  bkg_interval     numeric interval on which to alternate colors, 
#                           superceded by chk_bkg_idx
# @nparam  set_download     Have table include download link if set to a true
#                           value.  If text value is used then it will be used
#                           as text for the link.
# @nparam  sortable         If true, make data table sortable
#
#
# 
#
sub encodeSectionTable {
  my $METHOD = 'encodeSectionTable';
  my $self = shift || die ("self not passed");
  my %args = @_;

  $args{rs_params} ||= {};

  my $pre_text = '';
  my $class_def = '';
  my $id = $args{table_id} || '';

  if ( $args{sortable} ) {
    if ( !$self->{_included_sortable} ) {
      $pre_text = qq~
      <SCRIPT LANGUAGE=javascript SRC="$HTML_BASE_DIR/usr/javascript/sorttable.js"></SCRIPT>
      ~;
      $self->{_included_sortable}++;
    }
    $class_def = $args{class} || 'PA_sort_table';
  }

  my @table_attrs = ( 'border' => 0, );
  my $tr_info = $args{tr_info} || 'NOOP=1 ';
  my $tab = SBEAMS::Connection::DataTable->new( @table_attrs,
						__use_thead => 1,
                                                __tr_info => $tr_info,
						CLASS => $class_def,
                                                ID => $id  );
  my $num_cols = 0;


  my $rs_link = '';
  my $rs_name = '';
  my $file_prefix = $args{file_prefix} || 'mrm_';

  my @rs_data =  @{$args{rows}};

  if ( $args{set_download} ) {

    my $rs_headers = shift( @rs_data );
    $rs_headers = $args{rs_headings} if $args{rs_headings};
    
    $rs_name = $self->make_resultset( rs_data => \@rs_data, 
                                      headers => $rs_headers,
                                  file_prefix => $file_prefix,
                                   rs_params => $args{rs_params} );

    my $tsv_link = "<a href='$CGI_BASE_DIR/GetResultSet.cgi/$rs_name.tsv?rs_set_name=$rs_name&format=tsv;remove_markup=1' TITLE='Download table as tab-delimited text file'>TSV</a>";
    my @downloads = ( $tsv_link ); 
    if ( $args{download_links} ) {
      push @downloads, @{$args{download_links}};
    }
    my $download_links = join( ', ', @downloads );
    $rs_link = "<SPAN CLASS=info_box>Download as: $download_links</SPAN>";

    if ( $args{download_form} ) {
      my $hidden = qq~
      <INPUT TYPE=HIDDEN NAME=rs_set_name VALUE=$rs_name>
      <INPUT TYPE=HIDDEN NAME=remove_markup VALUE=1>
      <INPUT TYPE=HIDDEN NAME=tsv_output VALUE=1>
      ~;
      $rs_link = $args{download_form};
      $rs_link =~ s/HIDDEN_PLACEHOLDER/$hidden/m;
    }

  }

  return '' unless $args{rows};
  $args{header} ||= 0;
  if ( $args{width} ) {
    push @table_attrs, 'WIDTH', $args{width};
  }

  $args{max_rows} ||= 0;
  my $sbeams = $self->getSBEAMS();
  my $prefix = $sbeams->getRandomString( num_chars => 8, 
                                          char_set => ['A'..'Z', 'a'..'z'] );
  my $first = 1;
  my $bgcolor = '#f3f1e4'; #C0D0C0';
  my $chg_idx;
  my $rcnt = ( $args{header} ) ? 1 : 0;
#
  my $msg = '';
  for my $row ( @{$args{rows}} ) {
    $num_cols = scalar( @$row ) unless $num_cols;
    $tab->addRow( $row );
    $rcnt++;
    if ( defined $args{chg_bkg_idx} ) { # alternate on index
      if ( !$chg_idx ) {
        $chg_idx = $row->[$args{chg_bkg_idx}];
      } elsif ( $chg_idx ne $row->[$args{chg_bkg_idx}] ) {
#       $bgcolor = ( $bgcolor eq '#C0D0C0' ) ? '#F5F5F5' : '#C0D0C0';
        $bgcolor = ( $bgcolor eq '#f3f1e4' ) ? '#d3d1c4' : '#f3f1e4';
        $chg_idx = $row->[$args{chg_bkg_idx}];
      }

    } elsif ( $args{bkg_interval} ) { # alternate on n_rows
	unless ( $rcnt % $args{bkg_interval} ) {
#        $bgcolor = ( $bgcolor eq '#C0D0C0' ) ? '#F5F5F5' : '#C0D0C0';
	  $bgcolor = ( $bgcolor eq '#f3f1e4' ) ? '#d3d1c4' : '#f3f1e4';
      }
    } elsif ( $args{bg_color} ) { # single solid color
      $bgcolor = $args{bg_color};
    }
    $tab->setRowAttr(  ROWS => [$tab->getRowNum()], BGCOLOR => $bgcolor );
    my $num_rows = $tab->getRowNum() - 1;

    if ( $args{rows_to_show} && $args{rows_to_show} < $num_rows ) {
      $tab->setRowAttr( ROWS => [$tab->getRowNum()], ID => $prefix . '_toggle', 
                                                     NAME => $prefix . '_toggle', 
                                                     CLASS => 'hidden' ); 
    }
    # Message regarding truncated rows.
    if ( $args{max_rows} && $args{max_rows} <= $num_rows ) {
      my $span = scalar( @$row );
      $msg = "Table truncated at $args{max_rows} rows";
      if  ( $args{set_download} ) {
        $msg .= "(won't affect download, " . scalar( @{$args{rows}} ) . ' total rows)';
      }
      if ( $args{truncate_msg} ) {
        $msg .= ". $args{truncate_msg}";
      }
      if ( !$args{truncate_msg_as_text} ) {
        $tab->addRow( [$sbeams->makeInfoText("<I>$msg</I>")] );
      }
      $tab->setColAttr( ROWS => [$tab->getRowNum()], 
                     COLSPAN => $span, 
                        COLS => [ 1 ], 
                       ALIGN => 'CENTER' ); 

      if ( $args{rows_to_show} && $args{rows_to_show} < $num_rows ) {
        $tab->setRowAttr( ROWS => [$tab->getRowNum()], ID => $prefix . '_toggle', 
                                                     NAME => $prefix . '_toggle', 
                                                     CLASS => 'hidden' ); 
      } 
      last;
    }
  }

  my $nowrap = $args{nowrap} || [];
  if ( scalar( @$nowrap ) ) {
    $tab->setColAttr( COLS => $nowrap, NOWRAP => 1 ); 
  }
  
  # How many do we have?
  my $tot = $tab->getRowNum();
  my $closelink;
  if ( $args{rows_to_show} && $args{rows_to_show} < $tot - 1 ) {
    $closelink = $self->add_tabletoggle_js(); 
    $closelink .= "\n<FONT COLOR=BLUE><A HREF=#null ONCLICK=toggle_em('$prefix');return><SPAN ID='${prefix}_text' NAME='${prefix}_text' >[Show more rows]</A></FONT>";
  }

#  # if no wrapping desired...
#  if ( $args{nowrap} ) {
#    $tab->setRowAttr( ROWS => [1..$tot], NOWRAP => 1 ); 
#  }

  # Set header attributes
  if ( $args{header} ) {
    if ( $args{sortable} ) {
#      $tab->setRowAttr( ROWS => [1], BGCOLOR => '#0000A0', CLASS => 'sortheader' );
      $tab->setRowAttr( ROWS => [1], BGCOLOR => '#003F72', CLASS => 'sortheader' );
    } else {
      $tab->setRowAttr( ROWS => [1], BGCOLOR => '#003F72', CLASS => 'sortheader' );
    }
    $tab->setRowAttr( ROWS => [1], ALIGN => 'CENTER' ) unless $args{nocenter};
  }

  if ( $args{align} ) {
    for ( my $i = 0; $i <= $#{$args{align}}; $i++ ) {
      $tab->setColAttr( ROWS => [2..$tot], COLS => [$i + 1], ALIGN => $args{align}->[$i] );
    }
  }
#  $tab->addRow( [$closelink] );
  if ( $args{table_only} ) {
    return "$tab";
  }

  my $html = "$pre_text\n";
  my $help = $args{help_text} || '';

  if ( $html && $args{manual_widgets} ) {
    $html .= "<tr $tr_info><td nowrap ALIGN=left>$closelink</td><td align=center>$help</td><td nowrap ALIGN=right>$rs_link</td></tr>\n";
  }

  if ( 0 && ( $rs_link || $args{change_form} ) ) {
    
    if ( !$rs_link ) {
      $html .= "<TR><TD NOWRAP ALIGN=left>$args{change_form}</TD></TR>\n";
    } elsif ( !$args{change_form} ) {
      $html .= "<TR><TD NOWRAP COLSPAN=$num_cols ALIGN=right>$rs_link</TD></TR>\n";
    } else {
      $html .= "<TR><TD NOWRAP ALIGN=left>$args{change_form}</TD>\n";
      $html .= "<TD NOWRAP ALIGN=right>$rs_link</TD></TR>\n";
    }
  }

  my $colspan = $args{colspan} || 2;
  $html .= "<TR><TD NOWRAP COLSPAN=$colspan>$tab</TD></TR>";
  $html .= '</TABLE>' if $args{close_table};
   
  if ( wantarray ) {

    if ( $args{truncate_msg_as_text} ) {
      $html = "<tr $tr_info><td>$closelink</td></tr>" . $html;
      return ($html, $rs_name, $msg);
    } else {
      return ($html, $rs_name);
    }
  } else { 
    if ( $args{unified_widgets} ) {
      my $widget = "<tr $tr_info><td align=left>$closelink</td>";
      $widget .= "<td align=center>$help</td>";
      $widget .= "<td align=right nowrap=1>$rs_link</td>"; 
      $widget .= "</tr>";
      return( $widget . $html ); 
    } else {
      $html = "<tr $tr_info><td>$closelink</td></tr>" . $html;
      return $html;
    }
  }
}

###############################################################################
# display sample plot
###############################################################################
sub getSamplePlotDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'getSamplePlotDisplay';

  for my $arg ( qw( n_obs obs_per_million ) ) {
    unless ( defined ($args{$arg} ) ) {
      $log->error( "Missing required argument $arg" );
      return undef;
    }
  }

  $args{no_header} ||= 0;

  my $header = '';
  unless ( $args{no_header} ) {

    if ( $args{link} ) {
      $header = $self->encodeSectionHeader( text => 'Observed in Experiments:',
					  anchor => 'samples',
					    link => $args{link} , LMTABS=>1, no_toggle=>1);
    } else {
      $header = $self->encodeSectionHeader( text => 'Observed in Experiments:',
					  anchor => 'samples', LMTABS=>1, no_toggle=>1 );
    }
  }

  my $trinfo = $args{tr_info} || '';
  my $height = 50  + scalar( @{$args{n_obs}} ) * 12;

  # unable to hide legend as desired. Works on google viz playground, but not here (Alas!).
  my $GV = SBEAMS::Connection::GoogleVisualization->new();
  my $chart = $GV->setDrawBarChart(  samples => $args{n_obs},
                                     options => '', # qq~ legend: { position: 'none'}, title: "Number of Observations" ~, 
                                  data_types => [ 'string', 'number' ],
                                    headings => [ 'Sample Name', 'Number of Observations' ],
                                    chart_div => 'chart1_div',
                                  );
  my $chart_2 = $GV->setDrawBarChart(samples => $args{obs_per_million},
                                     options => '', # qq~ legend: { position: 'none'}, title: "Obs per million spectra" ~,
                                  data_types => [ 'string', 'number' ],
                                    headings => [ 'Sample Name', 'Obs per million spectra' ],
                                   chart_div => 'chart1_div',
                                      no_div => 1,

      );
  my $h_info = $GV->getHeaderInfo();
  $chart = qq~
    <script type='text/javascript'>
    function toggle_plot() { 
      if ( window['next_plot'] == undefined ) {
        window['next_plot'] = 'chart1_div';
      }
      if ( window['next_plot'] == 'chart1_div' ) {
        draw_chart1();
        window['next_plot'] = 'chart2_div';
        document.getElementById('toggle_button').innerHTML = 'Show Obs Per Million Spectra Plot';
      } else {
        draw_chart2();
        window['next_plot'] = 'chart1_div';
        document.getElementById('toggle_button').innerHTML = 'Show Total Obs Plot';
      }
    } 
  </script>
    $h_info
    <TR $trinfo>
      <TD></TD>
      <TD><button type='button' id='toggle_button' onclick=toggle_plot()>Show Total Obs Plot</button> 
          &nbsp;$chart_2
          &nbsp;$chart 
      </TD>
    </TR>
  ~;

  return ( wantarray() ) ? ($header, $chart) : $header . "\n" . $chart;
}

###############################################################################
# displaySampleMap
###############################################################################
sub getSampleMapDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = ( peptide_field => 'peptide_accession', @_ );

  my $in = join( ", ", keys( %{$args{instance_ids}} ) );
  return unless $in;

  my $header = '';
  if ( $args{link} ) {
    $header .= $self->encodeSectionHeader( text => 'Experiment peptide map:',
					   link => $args{link},
					   anchor => 'samples'
                                          );
  } else {
    $header .= $self->encodeSectionHeader( text => 'Experiment peptide map:',
					   anchor => 'samples'
                                         );
  }
  $header = '' if $args{no_header};

  my $html = '';
  my $trinfo = $args{tr_info} || '';


  my $sql = qq~     
  	SELECT DISTINCT SB.atlas_search_batch_id, sample_tag, 
		PISB.n_observations, $args{peptide_field},
    CASE when n_genome_locations = 1 THEN 1 ELSE 2 END
		FROM $TBAT_ATLAS_SEARCH_BATCH SB 
	  JOIN $TBAT_SAMPLE S ON s.sample_id = SB.sample_id
	  JOIN $TBAT_PEPTIDE_INSTANCE_SEARCH_BATCH PISB ON PISB.atlas_search_batch_id = SB.atlas_search_batch_id
	  JOIN $TBAT_PEPTIDE_INSTANCE PI ON PI.peptide_instance_id = PISB.peptide_instance_id
	  JOIN $TBAT_PEPTIDE P ON P.peptide_id = PI.peptide_id
    WHERE PI.peptide_instance_id IN ( $in )
    AND S.record_status != 'D'
    -- ORDER BY PISB.n_observations, $args{peptide_field} ASC
    ORDER BY sample_tag ASC
  ~;

  my @samples = $sbeams->selectSeveralColumns($sql);

  my $sample_js;
  my %samples;
  my %peptides;
  my $cntr = 0;
  for my $row ( @samples ) { 
    $cntr++;
    my $key = $row->[1] . '::::' . $row->[0];
    $row->[3] .= '*' if $row->[4] < 2;
    $peptides{$row->[3]}++;
    $samples{$key} ||= {};
    $samples{$key}->{$row->[3]} = $row->[2];
  }
  my $array_def = qq~
    <script type="text/javascript">
    google.setOnLoadCallback(drawHeatMap);
    function drawHeatMap() {
    var data = new google.visualization.DataTable();
    data.addColumn('string', 'Sample Name');
    ~;
  my @peps = ( $args{force_order} ) ? @{$args{force_order}} : sort( keys( %peptides ) );
#  for my $pa ( sort( keys( %peptides ) ) ) {
  for my $pa ( @peps ) {
    $array_def .= "    data.addColumn('number', '$pa');\n";
  }
  $array_def .= "DEFINE_ROWS_HERE\n";

  my $row = 0;
  my $max = 0;
  my $min = 50;
  $cntr = 0;
  for my $sa ( sort { "\L$a" cmp "\L$b" } ( keys( %samples ) ) ) {
    my $col = 0;
    my ( $name, $id ) = split "::::", $sa;
    $array_def .= "    data.setValue( $row, $col, '$name' );\n";
    $col++;

    for my $pa ( @peps ) {
      if ( $samples{$sa}->{$pa} ) {
        my $pep_cnt = log(1 + $samples{$sa}->{$pa})/log(10);
	$max = ( $pep_cnt > $max ) ? $pep_cnt : $max;
	$min = ( $pep_cnt < $min ) ? $pep_cnt : $min;
        $array_def .= "    data.setValue( $row, $col, $pep_cnt );\n";
	$cntr++;
      }
      $col++;
    }
    $row++;
  }
  $array_def =~ s/DEFINE_ROWS_HERE/data.addRows($row);/;
  my $num_colors = 256;
  $array_def .= qq~
    heatmap = new org.systemsbiology.visualization.BioHeatMap(document.getElementById('heatmapContainer'));
    heatmap.draw(data, {numberOfColors:$num_colors,passThroughBlack:false,startColor:{r:255,g:255,b:255,a:1},endColor:{r:100,g:100,b:100,a:1},emptyDataColor:{r:256,g:256,b:256,a:1}});
		}
  </script>
  ~;

  $args{header_text} = ( $args{header_text} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=section_description>$args{header_text}</TD></TR>" : '';
  $args{second_header} = ( $args{second_header} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=info_text>$args{second_header}</TD></TR>" : '';
  my $content = qq~
  <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/prototype/1.6.1.0/prototype.js"></script>
  <script type="text/javascript" src="https://www.google.com/jsapi"></script>
  <script type="text/javascript">
    google.load("visualization", "1", {});
  </script>
  <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/main/js/load.js"></script>
  <script type="text/javascript">
    systemsbiology.load("visualization", "1.0", {packages:["bioheatmap"]});
  </script>
  $array_def

  $args{header_text}
  $args{second_header}
	<TR $trinfo><TD> <DIV ID="heatmapContainer"></DIV>  </TD></TR>
	~;

  return ( wantarray() ) ? ($header, $content) : $header . "\n" . $content;

} # end getSampleMapDisplay

sub getBuildSelector {
  my $self = shift;
  my %args = @_;
  my $build_id = $args{atlas_build_id};
  my $sbeams = $self->getSBEAMS();

  my $accessible_builds = join( ',', $self->getAccessibleBuilds() );
  my $accessible_projects = join( ',', $sbeams->getAccessibleProjects() );

  # Get a hash of available atlas builds
  my $sql = qq~
  SELECT atlas_build_id, atlas_build_name
  FROM $TBAT_ATLAS_BUILD
  WHERE project_id IN ( $accessible_projects )
  AND record_status!='D'
  ORDER BY atlas_build_name
  ~;
  my @ordered_ids;
  my %id2build;
  my $sth = $sbeams->get_statement_handle( $sql );
  while ( my @row = $sth->fetchrow_array() ) {
    push @ordered_ids, $row[0];
    $id2build{$row[0]} = $row[1];
    $build_id ||= $row[0];
  }
  my $build_selector =  $q->popup_menu( -name => "atlas_build_id",
                                        -values => [ @ordered_ids ],
                                        -labels => \%id2build,
                                        -default => $build_id,
                                        -onChange => 'switchAtlasBuild()' );
  my $selector_widget = qq~
    <form name=build_form id=build_form>
    $build_selector
    </form>
    <script LANGUAGE="Javascript">
      function switchAtlasBuild() {
        document.build_form.submit();
      }
    </script>
  ~;
  return $selector_widget;
}

sub getSampleMapDisplayMod {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = ( peptide_field => 'peptide_accession', @_ );

  my @all_peps = ( @{$args{snp_support}}, @{$args{snp_original}} );
  my @speps = @{$args{snp_support}};
  my $in = join( ", ", @all_peps );
  return unless $in;

  my $header = '';
  if ( $args{link} ) {
    $header .= $self->encodeSectionHeader( text => 'Sample peptide map:',
                                          link => $args{link},
                                          anchor => 'samples'
                                          );
  } else {
    $header .= $self->encodeSectionHeader( text => 'Sample peptide map:',
                                          anchor => 'samples'
                                         );
  }
  $header = '' if $args{no_header};

  my $html = '';
  my $trinfo = $args{tr_info} || '';


  my $sql = qq~     
  	SELECT DISTINCT SB.atlas_search_batch_id, sample_tag, 
		PISB.n_observations, peptide_sequence, PI.peptide_instance_id
		FROM $TBAT_ATLAS_SEARCH_BATCH SB 
	  JOIN $TBAT_SAMPLE S ON s.sample_id = SB.sample_id
	  JOIN $TBAT_PEPTIDE_INSTANCE_SEARCH_BATCH PISB ON PISB.atlas_search_batch_id = SB.atlas_search_batch_id
	  JOIN $TBAT_PEPTIDE_INSTANCE PI ON PI.peptide_instance_id = PISB.peptide_instance_id
	  JOIN $TBAT_PEPTIDE P ON P.peptide_id = PI.peptide_id
    WHERE PI.peptide_instance_id IN ( $in )
    AND S.record_status != 'D'
    ORDER BY sample_tag ASC
  ~;

  my @samples = $sbeams->selectSeveralColumns($sql);

  my $sample_js;
	my %samples;
	my %peptides;
	my $cntr = 0;
  my %id2seq;
  for my $row ( @samples ) { 
		$cntr++;
		my $key = $row->[1] . '::::' . $row->[0];
		$peptides{$row->[3]}++;
		$samples{$key} ||= {};
		$samples{$key}->{$row->[3]} = $row->[2];
    $id2seq{$row->[4]} = $row->[3];
	}
	my $array_def = qq~
	<script type="text/javascript">
    google.setOnLoadCallback(drawHeatMap);
    function drawHeatMap() {
    var data = new google.visualization.DataTable();
    data.addColumn('string', 'Sample Name');
	~;
#	my @peps = ( $args{force_order} ) ? @{$args{force_order}} : sort( keys( %peptides ) );
  my %speps;
  for my $pep ( @speps ) {
    $speps{$pep}++;
  }

#  for my $pa ( sort( keys( %peptides ) ) ) {
  for my $pi ( @all_peps ) {
    my $pa = $id2seq{$pi};
    $array_def .= "    data.addColumn('number', '$id2seq{$pi}');\n";
	}
	$array_def .= "DEFINE_ROWS_HERE\n";

	my $row = 0;
	my $max = 0;
	my $min = 50;
	$cntr = 0;
	for my $sa ( sort( keys( %samples ) ) ) {
	  my $col = 0;
		my ( $name, $id ) = split "::::", $sa;
    $array_def .= "    data.setValue( $row, $col, '$name' );\n";
	  $col++;

    for my $pi ( @all_peps ) {
      my $pa = $id2seq{$pi};
			if ( $samples{$sa}->{$pa} ) {
        my $pep_cnt = log(1 + $samples{$sa}->{$pa})/log(10);
			  $max = ( $pep_cnt > $max ) ? $pep_cnt : $max;
    		$min = ( $pep_cnt < $min ) ? $pep_cnt : $min;
        $array_def .= "    data.setValue( $row, $col, $pep_cnt );\n";
		    $cntr++;
			}
		  $col++;
		}
		$row++;
	}
	$array_def =~ s/DEFINE_ROWS_HERE/data.addRows($row);/;
	my $num_colors = 256;
	$array_def .= qq~
	heatmap = new org.systemsbiology.visualization.BioHeatMap(document.getElementById('heatmapContainer'));
	heatmap.draw(data, {numberOfColors:$num_colors,passThroughBlack:false,startColor:{r:255,g:255,b:255,a:1},endColor:{r:100,g:100,b:100,a:1},emptyDataColor:{r:256,g:256,b:256,a:1}});
		}
  </script>
  ~;

	$args{header_text} = ( $args{header_text} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=section_description>$args{header_text}</TD></TR>" : '';
	$args{second_header} = ( $args{second_header} ) ? "<TR $trinfo><TD ALIGN=CENTER CLASS=info_text>$args{second_header}</TD></TR>" : '';
	my $content = qq~
  <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/prototype/1.6.1.0/prototype.js"></script>
  <script type="text/javascript" src="https://www.google.com/jsapi"></script>
  <script type="text/javascript">
    google.load("visualization", "1", {});
  </script>
  <script type="text/javascript" src="$HTML_BASE_DIR/usr/javascript/main/js/load.js"></script>
  <script type="text/javascript">
    systemsbiology.load("visualization", "1.0", {packages:["bioheatmap"]});
  </script>

	$array_def

  $args{header_text}
  $args{second_header}
	<TR $trinfo><TD> <DIV ID="heatmapContainer"></DIV>  </TD></TR>
	~;

  return ( wantarray() ) ? ($header, $content) : $header . "\n" . $content;


} # end getSampleMapDisplay



### END

###############################################################################
# display enhance sample info
###############################################################################
sub getPeptideSampleDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'getDetailedSampleDisplay';
  my $rows_to_show = $args{rows_to_show} || 25;
  my $mia = '';
  for my $arg ( qw ( sample_ids build_clause peptide_clause ) ) {
    next if defined $args{$arg};
    my $sep = ( $mia ) ? ',' : '';
    $mia .= $mia . $sep . $arg;
  }
  $args{tr_info} ||= '';
  if ( $mia ) {
    $log->error( "Missing required argument(s} $mia" );
    return;
  }
  my $in = join( ", ", @{$args{sample_ids}} );
  return unless $in;

  my $sql = qq~
  	SELECT S.sample_id, sample_title, PISB.n_observations,
           instrument_name, CASE WHEN ENZ.name IS NULL THEN 'Trypsin' ELSE ENZ.name END AS Enzyme,
           PUB.publication_name, PUB.abstract , PUB.uri, S.repository_identifiers
		FROM $TBAT_ATLAS_SEARCH_BATCH SB 
	  JOIN $TBAT_SAMPLE S ON s.sample_id = SB.sample_id
	  JOIN $TBAT_PEPTIDE_INSTANCE_SEARCH_BATCH PISB ON PISB.atlas_search_batch_id = SB.atlas_search_batch_id
	  JOIN $TBAT_PEPTIDE_INSTANCE PI ON PI.peptide_instance_id = PISB.peptide_instance_id
	  JOIN $TBAT_PEPTIDE P ON P.peptide_id = PI.peptide_id
	  LEFT JOIN $TBPR_INSTRUMENT I ON S.instrument_model_id = I.instrument_id 
	  LEFT JOIN $TBAT_PROTEASES ENZ ON ENZ.id = S.protease_id
    LEFT JOIN $TBAT_SAMPLE_PUBLICATION SP ON SP.sample_id = S.sample_id  
    LEFT JOIN $TBAT_PUBLICATION PUB ON PUB.publication_id = SP.publication_id  
    WHERE S.sample_id IN ( $in )
    $args{build_clause}
    $args{peptide_clause}
    AND S.record_status != 'D'
    ORDER BY sample_title
  ~;
  my @rows = $sbeams->selectSeveralColumns($sql);
  my $table = $self -> getSampleTableDisplay (data => \@rows,
                  rows_to_show => $rows_to_show,
                  type => 'Peptide');
  return $table;
}


sub getSampleTableDisplay{
  my $self = shift;
  my %args = @_;
  my $data = $args{data};
  my $type = $args{type};
  my $rows_to_show = $args{rows_to_show} || 25;
  my @samples = ();
  my $pre_id = '';
  foreach my $row (@$data){
     my ($id, $title, $nobs,$ins,$enzyme,$pub_name, $abstract, $link, $rp_id) = @$row;
     $id = "<a href='". "$CGI_BASE_DIR/PeptideAtlas/ManageTable.cgi?TABLE_NAME=AT_SAMPLE&sample_id=" . $id . "' target='_blank'>$id</a>";
     $rp_id = $self->get_dataset_url ($rp_id);
     if ($pre_id eq $id){next;$pre_id = $id;}  ## keep one publication only
     if ($abstract){
       $pub_name = $self->make_pa_tooltip( tip_text => "Abstract: $abstract", link_text => "<a href='$link'>$pub_name</a>" );
     }
     if ($type eq 'Peptide'){
       push @samples, [$id, $rp_id, $title, $nobs,$ins,$enzyme,$pub_name];
     }else{
       push @samples, [$id, $rp_id, $title,$ins,$enzyme,$pub_name];
     }
     $pre_id = $id;
  }
  my @align=();
  my @headings;
  my @cols = ();
  if ($type eq 'Peptide'){
    @cols = ('Experiment ID','Dataset','Experiment Name','NObs','Instrument','Enzyme','Publication');
    @align = qw(center left left center left left);
  }else{
    @cols = ('Experiment ID','Dataset','Experiment Name', 'Instrument','Enzyme','Publication');
    @align = qw(center left left center left left left);
  }
  foreach my $col (@cols){
    push @headings, $col, '';
  }
  my $headings = $self->make_sort_headings( headings => \@headings,  default => 'Experiment ID');
  unshift @samples, ($headings);
  my $table = $self->encodeSectionTable( header => 1, 
                                         tr_info => $args{tr_info},
                                         align  => [@align],
																				 bkg_interval => 3,
                                         nowrap => [qw(4 6)],
                                         rows_to_show => $rows_to_show,
                                         max_rows => $args{max_rows},
                                         sortable => 1,
                                         rows => \@samples );
  return $table;
}


sub getPTMTableDisplay {
  my $self = shift;
  my %args = @_;
  my $SUB_NAME = 'getPTMTableDisplay';
  my $cols = $args{cols};
  my $data = $args{data};
  my $self_build_id = $args{atlas_build_id}; 
  my $biosequence_name = $args{biosequence_name};
  my @rows;
  shift @$cols;
  pop @$cols;
  push @rows, [@$cols];
  foreach my $prot (keys %$data){
    foreach my $pos (sort {$a <=> $b} keys %{$data->{$prot}}){
      my @row = ();
      push @row, $pos+1; 
      #next if ($data->{$prot}{$pos}{nObs} ==0 );
      foreach my $col (@$cols){
        next if ($col =~ /offset/i);
        if ($col eq 'Residue' && ($data->{$prot}{$pos}{nObs} > 0            ||
				  $data->{$prot}{$pos}{InUniprot}  eq 'yes' ||
				  $data->{$prot}{$pos}{InNextProt} eq 'yes' )) {
          my $start_in_biosequence = $pos + 1;
          my $link = "$CGI_BASE_DIR/PeptideAtlas/GetPeptide?".
                     "atlas_build_id=$self_build_id&searchWithinThis=Peptide+Sequence&searchForThis=".
                     "$data->{$prot}{$pos}{peptide}&apply_action=QUERY"; 
          $data->{$prot}{$pos}{$col} = $self->make_pa_tooltip( tip_text => "Get peptide sequence covering this site",
                                                               link_text => "<a href='$link'>$data->{$prot}{$pos}{$col}</a>" );
        }
        elsif ($col eq 'InUniprot'  && $data->{$prot}{$pos}{$col} eq 'no') {
	  $data->{$prot}{$pos}{$col} = '-';
	}
        elsif ($col eq 'InNextProt' && $data->{$prot}{$pos}{$col} eq 'no') {
          $data->{$prot}{$pos}{$col} = '-';
	}

        push @row, $data->{$prot}{$pos}{$col};
      }
      push @rows , [@row];
    }
  }
  my @align = ();
  foreach my $i(0..15){
    push @align, 'center';
  }
  my $table = $self->encodeSectionTable( header => 1,
                                         align  => [@align],
					 bkg_interval => 3,
                                         rows_to_show => $args{rows_to_show},
                                         max_rows => $args{max_rows},
                                         rows => \@rows );
  return $table;
}


###############################################################################
# get_individual_spectra_display 
###############################################################################
sub get_individual_spectra_display {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'get_individual_spectra_display';
  my $bg_color = $args{bg_color} || '';
  my $column_titles_ref = $args{column_titles_ref};
  my $colnameidx_ref = $args{colnameidx_ref};
  my $hidden_cols_ref = $args{hidden_cols_ref};

  my $sortable = 1;
  if ( $args{sortable} ne ''){
    $sortable = $args{sortable};
  }


  unless( $args{resultset_ref} ) {
    $log->error( "No resultset_ref passed to display" );
    return;
  }
  my $resultset_ref = $args{resultset_ref};
  my @data = ();
  my @headings =();

  my $loop =1;

  foreach my $row (@{$resultset_ref->{data_ref}}){
    my @tmp =();
    foreach my $col (sort {$colnameidx_ref->{$a} <=> $colnameidx_ref->{$b}} keys %$colnameidx_ref){
       if (not defined $hidden_cols_ref->{$col}){
         push @tmp, $row->[$colnameidx_ref->{$col}];
         push @headings, $column_titles_ref->[$colnameidx_ref->{$col}] if ($loop==1);
       }
    }
    $loop++;
    push @data, \@tmp;
  }
  my $headings = $self->get_column_defs( labels => \@headings, plain_hash => 1 );
  unshift @data, ($self->make_sort_headings(headings => $headings));

  
    #push @{$resultset_ref->{column_list_ref}}, 'num_prot_mappings';
    #push @{$resultset_ref->{types_list_ref}}, 'int';

  my $align = [qw(left center left left center left center center center left left left)];

  my $html = $self->encodeSectionTable( header => 1,
					unified_widgets => 1,
					set_download => 1,
					colspan => 3,
					align  => $align,
					rows => \@data,
					rows_to_show => 10,
					max_rows => 200,
					nowrap => [1],
					bkg_interval => 3,
					rs_headings => \@$column_titles_ref, 
					bg_color => '#EAEAEA',
					sortable => 1,
					close_table => 1,
					table_id => 'individual_spectra',
					truncate_msg => '  <b>Use link in Modified Peptides section to filter results</b>.',
      );

  return "<TABLE>$html\n";


} # end 

###############################################################################
# displaySamples
###############################################################################
sub getProteinSampleDisplay {
  my $self = shift;
  my $sbeams = $self->getSBEAMS();
  my %args = @_;
  my $SUB_NAME = 'getProteinSampleDisplay';
  my $bg_color = $args{bg_color} || '';
  my $sortable = 1;
  if ( $args{sortable} ne ''){
    $sortable = $args{sortable};
  }
  my $rows_to_show = $args{rows_to_show} || 25;

  unless( $args{sample_ids} ) {
    $log->error( "No samples passed to display samples" );
    return;
  }

  my $in = join( ", ", @{$args{sample_ids}} );
  return unless $in;
  my $sql = qq~
    SELECT S.SAMPLE_ID,
           S.SAMPLE_TITLE, 
           '' ,
           INSTRUMENT_NAME, 
           CASE WHEN ENZ.name IS NULL THEN 'Trypsin' ELSE ENZ.name END AS Enzyme, 
           PUB.PUBLICATION_NAME, 
           PUB.ABSTRACT , 
           PUB.URI, 
           S.REPOSITORY_IDENTIFIERS
    FROM $TBAT_SAMPLE S
    LEFT JOIN $TBPR_INSTRUMENT I ON S.instrument_model_id = I.instrument_id
    LEFT JOIN $TBAT_PROTEASES ENZ ON ENZ.id = S.protease_id
    LEFT JOIN $TBAT_SAMPLE_PUBLICATION SP ON SP.SAMPLE_ID = S.SAMPLE_ID
    LEFT JOIN $TBAT_PUBLICATION PUB ON PUB.PUBLICATION_ID = SP.PUBLICATION_ID
    WHERE S.SAMPLE_ID IN ( $in )
    AND S.RECORD_STATUS != 'D'
    ORDER BY sample_tag ASC
  ~;
  my @rows = $sbeams->selectSeveralColumns($sql);
  my $table = $self -> getSampleTableDisplay(data => \@rows,
                               rows_to_show => $rows_to_show, 
                               type => 'Protein');
  return $table;
} # end getProteinSampleDisplay 

sub add_tabletoggle_js {
  my $self = shift;
  return '' if $self->{_added_ttoggle_js};
  $self->{_added_ttoggle_js}++;
  return <<"  END";
  <STYLE TYPE="text/css" media="screen">
    tr.visible { display: block-row; }
    tr.hidden { display: none; }
    td.visible { display: table-cell; }
    td.hidden { display: none; }
  </STYLE>
  <SCRIPT TYPE="text/javascript">
    function toggle_em(prefix) {
      
      // Grab page elements by their IDs
      var togglekey = prefix + '_toggle';
      var textkey = prefix   + '_text';

      var rows = document.getElementsByName(togglekey);
      var show = document.getElementById(textkey);

      var ttext = show.innerHTML;
      if ( ttext == '[Show more rows]' ) {
        show.innerHTML = "[Show fewer rows]";
      } else {
        show.innerHTML = "[Show more rows]";
      }
      
      for (var i=0; i < rows.length; i++) {
        if ( rows[i].className == 'hidden' ) {
           rows[i].className = 'visible';
        } else {
           rows[i].className = 'hidden';
        }
      }
    }
    </SCRIPT>

  END
  
}

#+
# Routine will set up PA tooltip colors and javascript.  
#
# @return   JS/CSS to be printed in the page at an appropriate juncture
#-
sub init_pa_tooltip {
  my $self = shift;
  
  # return nothing if tooltip has already been initialized
  return '' if $self->{_tooltip_init};
  
  $self->{_tooltip_init} = 1;
  return <<"  END";
  <SCRIPT DEFER="defer" src="$HTML_BASE_DIR/usr/javascript/TipWidget.js" TYPE="text/javascript"></SCRIPT>
  <STYLE>
  div#tooltipID { background-color:#F0F0F0;
                  border:2px 
                  solid #FF8C00; 
                  padding:4px; 
                  line-height:1.5; 
                  width:auto; 
                  font-family: Helvetica, Arial, sans-serif;
                  font-size:12px; 
                  font-weight: normal; 
                  position:absolute; 
                  visibility:hidden; left:0; top:0;
                }
  </STYLE>
  END
}


#+
# Routine to generate a 'tooltip' mouseover widget with 
#-
sub make_pa_tooltip {
  my $self = shift;
  my %args = @_;
  for my $arg ( qw( tip_text link_text ) ) {
    if ( !$args{$arg} ) {
      $log->error( "Missing required argument $arg" );
      return "";
    }
  }

  # Issue warning if tooltip js wasn't printed (well, fetched anyway).
  $self->init_pa_tooltip();
  $log->warn( "Failed to itialize tooltip!" ) if !$self->{_tooltip_init};

  my $class = $args{class} || 'pseudo_link';

  # Sanitize tooltip text
  $args{tip_text} =~ s/\'/\\\'/g;
  $args{tip_text} =~ s/\"/\\\'/g;
  $args{tip_text} =~ s/\r\n/ /g;
  $args{tip_text} =~ s/\n/ /g;
  
  return "<SPAN CLASS=$class onMouseover=\"showTooltip(event, '$args{tip_text}')\" onMouseout=\"hideTooltip()\">$args{link_text}</SPAN>";
}

sub formatMassMods { 
  my $self = shift;
  my $sequence = shift || return undef;
  $sequence =~ s/\[/<SPAN CLASS="aa_mod">\[/gm;
  if ($sequence =~ /\)/){
    $sequence =~ s/\)/\)<\/SPAN>/gm;
    $sequence =~ s/\](\w|$)/\]<\/SPAN>$1/gm;
  }else{
    $sequence =~ s/\]/\]<\/SPAN>/gm;
  }
  return $sequence;
}


sub get_table_help_section {
  my $self = shift;
  my %args = @_;
  $args{showtext} ||= 'show column descriptions';
  $args{hidetext} ||= 'hide column descriptions';
  $args{heading} ||= '';
  $args{description} ||= '';
  $args{footnote} ||= '';


  my $ecnt = 0;
  my $index = "<TABLE class=info_box>\n";
  for my $entry ( @{$args{entries}} ) {
    $ecnt++;
    $index .= $self->encodeSectionItem( %$entry, nowrap_description => 1 );
  }
  $index .= "</TABLE>\n";

  my $content =<<"  END";
  <BR>
  <span class=section_heading>$args{heading}</span> 
  <span class=description>$args{description}</span>
  $index
  <span class=description>$args{footnote}</span>
  END

  $sbeams = $self->getSBEAMS();
  my $section_toggle = $sbeams->make_toggle_section( content => $content,
                                                     sticky   => 0,
                                                     visible  => 0,
                                                     imglink  => 1,
                                                     showimg  => "/info_small.gif",
                                                     hideimg  => "/info_small.gif",
                                                     textlink => 1,
                                                     name     => $args{name},
                                                     showtext => $args{showtext},
                                                     hidetext => $args{hidetext},
                                          );

  return "$section_toggle<BR>";
}


sub get_table_help {
  my $self = shift;
  my %args = @_;
  my $column_titles_ref = $args{column_titles_ref};
  my $colnameidx_ref = $args{colnameidx_ref};
  my $hidden_cols_ref = $args{hidden_cols_ref};
  $args{heading} ||= '';
  $args{description} ||= '';
  $args{footnote} ||= '';
  $args{showtext} ||= 'show column descriptions';
  $args{hidetext} ||= 'hide column descriptions';


  my @headings = ();
  if ($colnameidx_ref && $hidden_cols_ref){ 
		foreach my $col (sort {$colnameidx_ref->{$a} <=> $colnameidx_ref->{$b}} keys %$colnameidx_ref){
			 if (not defined $hidden_cols_ref->{$col}){
				 push @headings, $column_titles_ref->[$colnameidx_ref->{$col}];
			 }
		}
  }else{
   @headings = @$column_titles_ref;
  }

  my $headings = $self->get_column_defs( labels => \@headings );

  my $index = "<TABLE class=info_box>\n";
  for my $entry ( @$headings ) {
    $index .= $self->encodeSectionItem( %$entry, nowrap_description => 1 );
  }
  $index .= "</TABLE>\n";

  my $content =<<"  END";
  <BR>
  <span class=section_heading>$args{heading}</span> 
  <span class=description>$args{description}</span>
  $index
  <span class=description>$args{footnote}</span>
  END

  $sbeams = $self->getSBEAMS();
  my $section_toggle = $sbeams->make_toggle_section( content => $content,
                                                     sticky   => 0,
                                                     visible  => 0,
                                                     imglink  => 1,
                                                     showimg  => "/info_small.gif",
                                                     hideimg  => "/info_small.gif",
                                                     textlink => 1,
                                                     name     => $args{name},
                                                     showtext => $args{showtext},
                                                     hidetext => $args{hidetext},
                                          );

  return "$section_toggle<BR>";
}

###############################################################################
########### vocabHTML
######################################################################################

sub vocabHTML {

###### This subroutine displays the Vocab section on the HTML Page
 
my $self=shift;
my $sbeamsMOD=$self;

my ($table3);

my $pi_desc='Isoelectric point of the peptide';
my $ssrcalc_desc='Sequence Specific Retention Factor provides a hydrophobicity measure for each peptide using the algorithm of Krohkin et al. Version 3.0 <A HREF=http://hs2.proteome.ca/SSRCalc/SSRCalc.html target=_blank>more</A>';
my $org_desc='Organism in which this peptide is observed';
my $n_obs_desc='Number of MS/MS spectra that are identified to this peptide in each build. The hyperlink will take you to the peptide page that will display all the relevant information about this peptide contained in the listed PeptideAtlas build';
my $build_names_desc='Public build in which this peptide is found. The hyperlink will take you to a build specific page summarizing all the relevant information about the build.';

$table3 ='<BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><TABLE WIDTH=600>';

$table3 .= $sbeamsMOD->encodeSectionHeader( text => 'Vocabulary',
                                            width => 900
                                            );

#Substituting the path generated from encodeSectionHeader to match with peptideatlas.org for displaying orange header
$table3 =~ s/\/devAP\/sbeams//gm;

$table3 .= $sbeamsMOD->encodeSectionItem( key   => 'pI',
                                          value => $pi_desc,
                                          key_width => '20%'
                                         );
$table3.= $sbeamsMOD->encodeSectionItem( key   => 'SSRCalc',
                                          value =>$ssrcalc_desc
                                         );
$table3.= $sbeamsMOD->encodeSectionItem( key   => 'Organism Name',
                                          value =>$org_desc
                                         );

$table3.= $sbeamsMOD->encodeSectionItem( key   => 'Number of observations ',
                                          value =>$n_obs_desc
                                         );

$table3.= $sbeamsMOD->encodeSectionItem( key   => 'Build Names in which Peptide is found',
                                          value =>$build_names_desc
                                         );
$table3 .= '</TABLE>';

return ($table3);


}





sub get_atlas_checklist {
  my $self = shift;
  my %args = @_;

  # check input
	for my $arg ( qw( build_info js_call ) ) {
		return "Missing required arguement $arg" unless defined $args{$arg};
	}

  #### Get a list of accessible builds, use to validate passed list.
  my @accessible = $self->getAccessibleBuilds();
	my %build_check;
	for my $build ( @accessible ) {
		$build_check{$build}++;
	}

  # Convenience
	my %build_info = %{$args{build_info}};
  my @builds =  sort { $build_info{$a}->{org} cmp $build_info{$b}->{org} ||
		               $build_info{$b}->{is_curr} <=>  $build_info{$a}->{is_curr} ||
		            $build_info{$b}->{is_default} <=>  $build_info{$a}->{is_default} ||
						 		      $build_info{$a}->{name} cmp $build_info{$b}->{name} } keys %build_info;

  my $table = SBEAMS::Connection::DataTable->new(BORDER => 0);
  $table->addRow( [ 'Add', 'Build Id', 'Display Name', 'Full Name', 'Organism', 'is_def' ] );
  $table->setRowAttr(  ROWS => [1], BGCOLOR => '#bbbbbb', ALIGN=>'CENTER' );
  $table->setHeaderAttr( BOLD => 1 );

  for my $build ( @builds ) {
		my %build = %{$build_info{$build}};
    my $checked = ( $build{visible} ) ? 'checked' : '';

    my $chkbox =<<"    END";
    <INPUT $checked TYPE="checkbox" NAME="build_id" VALUE="$build" onchange="$args{js_call}($build);">
    END
    $table->addRow( [ $chkbox, $build, @build{qw(display_name name org is_default)} ] );
    $table->setRowAttr( ROWS => [$table->getRowNum()], BGCOLOR => $build{bgcolor} );

	}
	return "$table";
}

sub get_atlas_select {

  my $self = shift;
  my %args = @_;
  my $select_name = $args{select_name} || 'build_id';

  #### Get a list of accessible builds, use to validate passed list.
  my @accessible = $self->getAccessibleBuilds();
  my $build_string = join( ",", @accessible );
  return '' unless $build_string;

  my $sql = qq~
  SELECT atlas_build_id, build_tag 
  FROM $TBAT_ATLAS_BUILD
  WHERE atlas_build_id IN ( $build_string )
  ~;

  $log->info( $sql );

  my $select = "<INPUT TYPE=SELECT NAME=$select_name>\n";
  my $sth = $sbeams->get_statement_handle->( $sql );
  while ( my @row = $sth->fetchrow_array() ) {
    $select .= "<OPTION VALUE=$row[0]>$row[1]\n";
  }
  $select .= "</INPUT>\n";
  return $select;
}

sub get_proteome_coverage_new {
  my $self = shift;
  my $build_id = shift; 
  my $patterns = shift;
  my @patterns = @$patterns;

  return [] if (! $build_id || ! @patterns);
  my $sbeams = $self->getSBEAMS();
  my $sql = '';
  my $obs_sql = '';
  my @names = ();
  my $union = '';
  foreach my $line (@patterns){
    my ($org_id, $name, $type, $pat_str)  = split(/,/, $line);
    my @pats = split(/;/, $pat_str);
    my $contraint = '';
    my $or = '';
    push @names, $name;
    if ($type =~ /accession/i){
       foreach my $pat (@pats){
         $contraint .= "$or B2.BIOSEQUENCE_NAME LIKE '$pat%' ";
         $or = 'OR';
		   } 

    }elsif($type =~ /description/i){
       foreach my $pat (@pats){
         $contraint .= "$or B2.BIOSEQUENCE_Desc LIKE '%$pat%' ";
         $or = 'OR';
       }
     }
		 $sql .= qq~
			$union
			 (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
				 '$name' AS Name, 
				 B2.BIOSEQUENCE_SET_ID AS SETID
				FROM $TBAT_ATLAS_BUILD AB2
				JOIN $TBAT_BIOSEQUENCE B2 ON B2.BIOSEQUENCE_SET_ID = AB2.BIOSEQUENCE_SET_ID
				WHERE atlas_build_id = $build_id
				AND ($contraint) 
				GROUP BY B2.BIOSEQUENCE_SET_ID 
			 )
		 ~;
		$obs_sql .= qq~
			 $union
			(SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
						 '$name' AS CAT,
						 B2.BIOSEQUENCE_SET_ID AS SETID
				FROM $TBAT_BIOSEQUENCE B2
				JOIN $TBAT_BIOSEQUENCE_ID_ATLAS_BUILD_SEARCH_BATCH A
					ON A.BIOSEQUENCE_ID = B2.BIOSEQUENCE_ID
				JOIN $TBAT_ATLAS_BUILD_SEARCH_BATCH ABSB
					ON (A.ATLAS_BUILD_SEARCH_BATCH_ID = ABSB.ATLAS_BUILD_SEARCH_BATCH_ID)
				WHERE ABSB.atlas_build_id = $build_id
				AND ($contraint)
				GROUP BY B2.BIOSEQUENCE_SET_ID

			)
		 ~;
    $union = 'UNION';
  }

  $sql = qq~
  select * 
  FROM 
  (
    $sql
  ) As A 
  ~;
  #$log->error($sql);

  my $sth = $sbeams->get_statement_handle( $sql );
  my %entry_cnt = ();
  while ( my @row = $sth->fetchrow_array() ) {
    $entry_cnt{$row[1]} = $row[0];
  }
  my $obssth = $sbeams->get_statement_handle( $obs_sql );
  my %obs;
  while ( my @row = $obssth->fetchrow_array() ) {
    $obs{$row[1]} = $row[0];
  }
  my @headings;
  my %head_defs = ( Database => 'Name of database, which collectively form the reference database for this build',
                    N_Prots => 'Total number of entries in subject database',
                    N_Obs_Prots => 'Number of proteins within the subject database to which at least one observed peptide maps',
                    Pct_Obs => 'The percentage of the subject proteome covered by one or more observed peptides' );

  for my $head ( qw( Database N_Prots N_Obs_Prots Pct_Obs ) ) {
    push @headings, $head, $head_defs{$head};
  }
  my $headings = $self->make_sort_headings( headings => \@headings, default => 'Database' );
  my @return = ( $headings );

  for my $name ( @names ) {
    my $db = $name;
    my $obs  = $obs{$name} || 0;
    my $n_entry = $entry_cnt{$name} || 0;
    my $pct =0;
    if ( $obs && $n_entry ) {
        $pct = sprintf( "%0.1f", 100*($obs/$n_entry) );
    }
    push @return, [ $db, $n_entry, $obs, $pct ];
  }
  return '' if ( @return == 1);
  my $table = '<table width=600>';

  $table .= $self->encodeSectionHeader(
      text => 'Proteome Coverage (exhaustive)',
      no_toggle => 1,
      LMTABS => 1,
      width => 600
  );

  $table .= $self->encodeSectionTable( rows => \@return, 
                                        header => 1, 
                                        table_id => 'proteome_cover',
                                        align => [ qw(left right right right ) ], 
                                        bg_color => '#EAEAEA',
                                        rows_to_show => 25,
                                        sortable => 1 );
  $table .= '</TABLE>';

  return $table;
}

sub get_proteome_coverage {
  my $self = shift;
  my $build_id = shift || return [];

  my $sbeams = $self->getSBEAMS();

  my $sql = qq~
  select * 
  FROM(
		(
		SELECT COUNT(*) AS CNT, B.DBXREF_ID AS ID , DBXREF_NAME AS NAME , B.BIOSEQUENCE_SET_ID AS SETID
		FROM $TBAT_ATLAS_BUILD AB
		JOIN $TBAT_BIOSEQUENCE B ON B.biosequence_set_id = AB.biosequence_set_id
		JOIN biolink.dbo.dbxref DX ON DX.dbxref_id = B.dbxref_id
		WHERE atlas_build_id = $build_id
		AND B.dbxref_id IS NOT NULL
		GROUP BY B.dbxref_id, dbxref_name, B.biosequence_set_id
		)
		UNION 
		(SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT, 
					 '' AS ID, 
					 LEFT (B2.BIOSEQUENCE_NAME, PATINDEX('%[_]%', B2.BIOSEQUENCE_NAME)) AS NAME, 
					 B2.BIOSEQUENCE_SET_ID AS SETID
		FROM $TBAT_ATLAS_BUILD AB2
		JOIN $TBAT_BIOSEQUENCE B2 ON B2.BIOSEQUENCE_SET_ID = AB2.BIOSEQUENCE_SET_ID
		WHERE atlas_build_id = $build_id 
		AND B2.DBXREF_ID IS NULL
		AND B2.BIOSEQUENCE_NAME LIKE '%\\_%' ESCAPE '\\'
		AND B2.BIOSEQUENCE_NAME NOT LIKE 'DECOY%'
		AND B2.BIOSEQUENCE_NAME NOT LIKE 'ENST%'
     AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB%'
		GROUP BY B2.BIOSEQUENCE_SET_ID , LEFT (B2.BIOSEQUENCE_NAME, PATINDEX( '%[_]%', B2.BIOSEQUENCE_NAME ))
		)
    UNION
    (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
           '' AS ID,
           LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ),
           B2.BIOSEQUENCE_SET_ID AS SETID
    FROM $TBAT_ATLAS_BUILD AB2
    JOIN $TBAT_BIOSEQUENCE B2 ON B2.BIOSEQUENCE_SET_ID = AB2.BIOSEQUENCE_SET_ID
    WHERE atlas_build_id = $build_id
    AND B2.DBXREF_ID IS NULL
     AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB%'
    AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB_ENST%'
    GROUP BY B2.BIOSEQUENCE_SET_ID ,LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ) 
    )
    UNION
    (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
           '' AS ID,
           'CONTRIB_ENST_' AS NAME, 
           B2.BIOSEQUENCE_SET_ID AS SETID
    FROM $TBAT_ATLAS_BUILD AB2
    JOIN $TBAT_BIOSEQUENCE B2 ON B2.BIOSEQUENCE_SET_ID = AB2.BIOSEQUENCE_SET_ID
    WHERE atlas_build_id = $build_id
    AND B2.DBXREF_ID IS NULL
    AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB_ENST%'
    GROUP BY B2.BIOSEQUENCE_SET_ID
    )
  )CAT
  ORDER BY NAME
  ~;
  my $sth = $sbeams->get_statement_handle( $sql );
  my @names;
  while ( my @row = $sth->fetchrow_array() ) {
    push @names, \@row;
  }

  my $obs_sql = qq~
  (SELECT COUNT(DISTINCT B.BIOSEQUENCE_ID),convert(varchar(10), B.DBXREF_ID)
  FROM $TBAT_BIOSEQUENCE B
  JOIN $TBAT_PEPTIDE_MAPPING PM 
    ON PM.MATCHED_BIOSEQUENCE_ID = B.BIOSEQUENCE_ID
  JOIN $TBAT_PEPTIDE_INSTANCE PI 
    ON PM.PEPTIDE_INSTANCE_ID = PI.PEPTIDE_INSTANCE_ID
  WHERE atlas_build_id = $build_id
  GROUP BY B.DBXREF_ID
  )
  UNION
  (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT, 
         LEFT (B2.BIOSEQUENCE_NAME, PATINDEX('%[_]%', B2.BIOSEQUENCE_NAME)) AS CAT
    FROM $TBAT_BIOSEQUENCE B2
    JOIN $TBAT_PEPTIDE_MAPPING PM2
      ON PM2.MATCHED_BIOSEQUENCE_ID = B2.BIOSEQUENCE_ID
    JOIN $TBAT_PEPTIDE_INSTANCE PI2
      ON PM2.PEPTIDE_INSTANCE_ID = PI2.PEPTIDE_INSTANCE_ID
    WHERE atlas_build_id = $build_id
  and B2.dbxref_id is null
  AND B2.BIOSEQUENCE_NAME LIKE '%\\_%' ESCAPE '\\'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'DECOY%'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB%'
  GROUP BY LEFT (B2.BIOSEQUENCE_NAME, PATINDEX( '%[_]%', B2.BIOSEQUENCE_NAME ))
  )
  UNION
  (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
         LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ) AS CAT
    FROM $TBAT_BIOSEQUENCE B2
    JOIN $TBAT_PEPTIDE_MAPPING PM2
      ON PM2.MATCHED_BIOSEQUENCE_ID = B2.BIOSEQUENCE_ID
    JOIN $TBAT_PEPTIDE_INSTANCE PI2
      ON PM2.PEPTIDE_INSTANCE_ID = PI2.PEPTIDE_INSTANCE_ID
    WHERE atlas_build_id = $build_id
  and B2.dbxref_id is null
  AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB%'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB_ENST%'
  AND B2.BIOSEQUENCE_NAME NOT LIKE 'CONTRIB_H%'
  GROUP BY LEFT (B2.BIOSEQUENCE_NAME, charindex( '_' , B2.BIOSEQUENCE_NAME,charindex ('_',B2.BIOSEQUENCE_NAME, 0) +1) ) 
  )
  UNION
  (SELECT COUNT(DISTINCT B2.BIOSEQUENCE_ID) AS CNT,
         'CONTRIB_ENST_' AS CAT 
    FROM $TBAT_BIOSEQUENCE B2
    JOIN $TBAT_PEPTIDE_MAPPING PM2
      ON PM2.MATCHED_BIOSEQUENCE_ID = B2.BIOSEQUENCE_ID
    JOIN $TBAT_PEPTIDE_INSTANCE PI2
      ON PM2.PEPTIDE_INSTANCE_ID = PI2.PEPTIDE_INSTANCE_ID
    WHERE atlas_build_id = $build_id
  and B2.dbxref_id is null
  AND B2.BIOSEQUENCE_NAME LIKE 'CONTRIB_ENST%'
  )

  ~;

  my $obssth = $sbeams->get_statement_handle( $obs_sql );
  my %obs;
  while ( my @row = $obssth->fetchrow_array() ) {
    $obs{$row[1]} = $row[0];
  }

  my @headings;
  my %head_defs = ( Database => 'Name of databse, which collectively form the reference database for this build',
                    N_Prots => 'Total number of entries in subject database',
                    N_Obs_Prots => 'Number of proteins within the subject database to which at least one observed peptide maps',
                    Pct_Obs => 'The percentage of the subject proteome covered by one or more observed peptides' );

  for my $head ( qw( Database N_Prots N_Obs_Prots Pct_Obs ) ) {
    push @headings, $head, $head_defs{$head};
  }
  my $headings = $self->make_sort_headings( headings => \@headings, default => 'Database' );
  my @return = ( $headings );

  for my $row ( @names ) {
    my $db = $row->[2]; 
    $db =~ s/\_$//;
    if ( $db eq 'Swiss-Prot' ) {
      $db .= ' (may include Varsplic)';
    } elsif ( $db eq 'UniProt' ) {
      $db .= ' (excludes SwissProt if shown)';
    }
    my $pct = 0;
    my $obs = '';
    if($row->[1] =~ /^\d+$/ && $row->[1] !~ /^0$/){ 
      $obs = $obs{$row->[1]};
    }else{
      $obs = $obs{$row->[2]};
    }

    if ( $obs && $row->[0] ) {
        $pct = sprintf( "%0.1f", 100*($obs/$row->[0]) );
    }
    push @return, [ $db, $row->[0], $obs, $pct ];
  }
  return '' if ( @return == 1);
  my $table = '<table width=600>';

  $table .= $self->encodeSectionHeader(
      text => 'Proteome Coverage (exhaustive)',
      no_toggle => 1,
      LMTABS => 1,
      width => 600
  );

  $table .= $self->encodeSectionTable( rows => \@return, 
                                        header => 1, 
                                        table_id => 'proteome_cover',
                                        align => [ qw(left right right right ) ], 
                                        bg_color => '#EAEAEA',
                                        rows_to_show => 25,
                                        sortable => 1 );
  $table .= '</TABLE>';

  return $table;
}

sub get_what_is_new {
  my $self = shift;
  my $build_id = shift || return [];
  my $sbeams = $self->getSBEAMS();
  # check if it default build 
  my $sql = qq~
    SELECT  DEFAULT_ATLAS_BUILD_ID
    FROM $TBAT_DEFAULT_ATLAS_BUILD 
    WHERE ATLAS_BUILD_ID = $build_id
  ~;
  my @row = $sbeams->selectOneColumn($sql); 
  return if(! @row);
    
  # compare biosequence set version
  $sql = qq~
    SELECT AB.ATLAS_BUILD_NAME, BS.SET_DESCRIPTION    
    FROM $TBAT_ATLAS_BUILD AB
    JOIN $TBAT_BIOSEQUENCE_SET BS ON (BS.BIOSEQUENCE_SET_ID = AB.BIOSEQUENCE_SET_ID ) 
    WHERE AB.ATLAS_BUILD_ID = $build_id     
  ~;
  @row = $sbeams->selectSeveralColumns($sql);
  my ($default_build_name, $default_bsset_desc) = @{$row[0]};
  my $build_name_pat = $default_build_name;
  $build_name_pat =~ s/\s+\d.*//;
  $sql = qq~
    SELECT AB.atlas_build_id, AB.ATLAS_BUILD_NAME, BS.SET_DESCRIPTION
    FROM $TBAT_ATLAS_BUILD AB
    JOIN $TBAT_BIOSEQUENCE_SET BS ON (BS.BIOSEQUENCE_SET_ID = AB.BIOSEQUENCE_SET_ID )
    WHERE AB.ATLAS_BUILD_NAME like '$build_name_pat [0-9]%' and AB.project_id = 475 
    AND AB.ATLAS_BUILD_ID != $build_id
  ~;
  @row = $sbeams->selectSeveralColumns($sql);
  return if (! @row);
  my ($previous_build_id, $previous_build_name, $previous_bsset_desc) = @{$row[0]};

  my @headings;
  push @headings, '', '';
  push @headings, $default_build_name, '';
  push @headings ,  $previous_build_name, '';

  my $headings = $self->make_sort_headings( headings => \@headings, default => $default_build_name );
  my @return = ( $headings );
  push @return , ['Reference_Database', $default_bsset_desc, $previous_bsset_desc];

  ## get sample number
  $sql = qq~ 
      SELECT ATLAS_BUILD_ID , COUNT(SAMPLE_ID)  
      FROM $TBAT_ATLAS_BUILD_SAMPLE 
      WHERE atlas_build_id in ($build_id, $previous_build_id)
      GROUP BY ATLAS_BUILD_ID
  ~;
  my %sample_cnt = $sbeams->selectTwoColumnHash($sql);


  $sql = qq~
  SELECT  atlas_build_id,  COUNT(distinct PI.peptide_instance_id) cnt
  FROM $TBAT_PEPTIDE_INSTANCE PI
  JOIN $TBAT_PEPTIDE_MAPPING PM ON (PI.PEPTIDE_INSTANCE_ID = PM.PEPTIDE_INSTANCE_ID)
  JOIN $TBAT_BIOSEQUENCE B ON (PM.MATCHED_BIOSEQUENCE_ID = B.BIOSEQUENCE_ID)
  WHERE atlas_build_id in ($build_id, $previous_build_id)
  AND B.BIOSEQUENCE_NAME NOT LIKE 'DECOY%'
  AND B.BIOSEQUENCE_NAME NOT LIKE 'CONTAM%'
  AND B.biosequence_name NOT LIKE '%UNMAPPED%'
  AND B.biosequence_desc NOT LIKE '%common contaminant%'
  GROUP BY atlas_build_id
  ~;
  my %pep_count = $sbeams->selectTwoColumnHash($sql);

  $sql = qq~ 
  SELECT PID.atlas_build_id, COUNT(BS.biosequence_name) cnt
  FROM $TBAT_PROTEIN_IDENTIFICATION PID
  JOIN $TBAT_PROTEIN_PRESENCE_LEVEL PPL
  ON PPL.protein_presence_level_id = PID.presence_level_id
  JOIN $TBAT_BIOSEQUENCE BS
  ON BS.biosequence_id = PID.biosequence_id
  WHERE PID.atlas_build_id in ($build_id, $previous_build_id)
  AND PPL.level_name = 'canonical'
  AND BS.biosequence_name NOT LIKE 'DECOY%'
  AND BS.biosequence_name NOT LIKE '%UNMAPPED%'
  AND BS.biosequence_desc NOT LIKE '%common contaminant%'
  AND BS.BIOSEQUENCE_NAME NOT LIKE 'CONTAM%'
  GROUP BY PID.atlas_build_id
  ~;
  my %prot_count = $sbeams->selectTwoColumnHash($sql);

  push @return , ['# Experiments', $sample_cnt{$build_id}, $sample_cnt{$previous_build_id}];
  push @return , ['Distinct_Peptides', $pep_count{$build_id}, $pep_count{$previous_build_id}];
  push @return , ['Canonical_Proteins', $prot_count{$build_id}, $prot_count{$previous_build_id}];

  my $table = '<table width="600">';
  $table .= $self->encodeSectionHeader(
      LMTABS => 1,
      no_toggle => 1,
      text => "What&#39s new",
  );

  $table .= $self->encodeSectionTable( rows => \@return,
                                        header => 1,
                                        table_id => 'what_is_new',
                                        align => [ qw(left right right right ) ],
                                        bg_color => '#EAEAEA',
                                        rows_to_show => 25,
                                        sortable => 1 );
  $table .= '</table>';
  ## new sample table:
  $sql = qq~
    SELECT SAMPLE_ID
    FROM $TBAT_SAMPLE 
    WHERE SAMPLE_ID 
    IN ( 
      SELECT A.SAMPLE_ID FROM $TBAT_ATLAS_BUILD_SAMPLE A 
      WHERE A.ATLAS_BUILD_ID IN ($build_id)
    ) 
    AND SAMPLE_ID 
    NOT IN ( 
      SELECT B.SAMPLE_ID 
      FROM $TBAT_ATLAS_BUILD_SAMPLE B 
      WHERE B.ATLAS_BUILD_ID IN ($previous_build_id)
    ) 
  ~;
  my @sample_ids = $sbeams->selectOneColumn($sql);
  if(@sample_ids){
    my ( $tr, $link ) = $sbeams->make_table_toggle( name => 'getprotein_samplemap',
                                                  visible => 1,
                                                tooltip => 'Show/Hide Section',
                                                 imglink => 1,
                                                  sticky => 1 );

    my $sampleDisplay = $self->getProteinSampleDisplay( sample_ids => \@sample_ids,
                                                          'link' => $link,
                                                     rows_to_show => 5,
                                                         max_rows => 500,
                                                         bg_color  => '#EAEAEA',
                                                         sortable => 1,
                                                         tr_info => $tr );

    $table .= "<TABLE width='600'>$sampleDisplay</TABLE>";
  }
  
  return $table;

}

sub get_scroll_table {
  my $self = shift || die ("self not passed");
  my %args = ( width => 900,
               bold => 1,
               key_width => 20,
               @_ );

  for my $arg ( qw( sql headings ) ) {
    die "Missing required argument $arg" if !defined $args{$arg};
  }

  $sbeams = $self->getSBEAMS();
  my $sth = $sbeams->get_statement_handle->( $args{sql} );
  while ( my @row = $sth->fetchrow_array() ) {
  }

}

sub print_html_table_to_tsv{
  my $self = shift;
  my %args = @_;
  my $data_ref = $args{data_ref};
  my $column_name_ref = $args{column_name_ref};
  my $filename = $args{filename};
  my %params =();
  $params{'Content-Disposition'}="attachment;filename=$args{filename}";
  print $q->header(-type=>'tsv',%params);

  foreach (@$column_name_ref){
    print lc($_) ."\t";
  }
  print "\n";
  foreach my $row (@$data_ref){
    foreach my $val (@$row){
      $val =~ s/(\<[^\<]+\>)/ /g;
      print "$val\t";
    }
    print "\n";
  }

}
sub display_peptide_sample_category_plotly{
  my $self = shift;
  my %args = @_;
  my $data_ref = $args{data_ref};
  my $sample_arrayref = $args{sample_arrayref};
  my $build_id = $args{build_id};
  my $column_name_ref = $args{column_name_ref};

  my %cols=();
  my $idx = 0;
  foreach (@$column_name_ref){
    $cols{$_} = $idx;
    $idx++;
  }
  my (@sample_category, @peptide_count,@links, @obs_per_million );
  my $total_observed_spectra;
  my @n_good_spectra; 

 foreach my $data (@$data_ref){
    my ($name,$id,$cnt) = @$data;
    my ($good_spectra) =0;
    foreach my $row (@$sample_arrayref){
      my $n = scalar @$row; 
      my $sample_cat_id  = $row->[$cols{"sample_category_id"}];
      if ($sample_cat_id eq $id){
        $good_spectra +=  $row->[$cols{"Spectra ID'd"}];
        $total_observed_spectra += $row->[$cols{"Spectra ID'd"}]; 
      }
    }  
 
    push @n_good_spectra, $good_spectra;
    push @peptide_count, $cnt; 
    push @sample_category,'<a href="'."https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptides?atlas_build_id=$build_id&sample_category_id=$id&QUERY_NAME=AT_GetPeptides&apply_action=QUERY" .'">'. $name .'  </a>';
  }

  foreach my $val (@peptide_count){
    push @obs_per_million,  sprintf( "%0.1f",($val/$total_observed_spectra) * 1000000 );
  }
  
  my $sample_category_str = join("','", @sample_category);
  my $n_good_spectra_str = join(",", @n_good_spectra);  
  my $obs_per_million_str = join(",", @obs_per_million);
  my $nc = $#sample_category + 1;
  my $height = $nc*40 > 1100? 1100:$nc*40; 
  my $plot_js = qq~
			var pepcnt = {
				y: ['$sample_category_str'],
				x: [$obs_per_million_str],
        //customdata:['link_str'],
				name:'Number Distinct Peptide Per Million Observed Spectra',
			  type: 'bar',
        marker: {color: '0x5588bb'},	
        orientation: 'h',
        opacity: 0.8
			};
      var obs = {
        y: ['$sample_category_str'],
        x: [$n_good_spectra_str],
        name:'# Spectrua ID',
        type: 'bar',
        marker: {color: '0x5588bb'},
        orientation: 'h',
        opacity: 0.8,
      };

			var layout = {
        legend:{x: 0.029,y: 1.1},
			  height: $height,
        font: {size: 16},
				title:'Number of Distinct Peptides Per Million Observed Spectra',
        margin:{l: 350},
        hoverlabel:{bgcolor:'white'},
			};
			var data1 = [pepcnt];
			Plotly.newPlot('plot_div3', data1,layout);
      layout.title = "Total Observed Spectra";
      var data2 = [obs];
      Plotly.newPlot('plot_div4', data2,layout);

  ~;
  my $chart = qq~
    <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
	  <script type='text/javascript'>
			function toggle_plot() {
        var myplot3 =  document.getElementById('plot_div3');
        var myplot4 =  document.getElementById('plot_div4');
        if (myplot3.style.display == '' ){
           myplot3.style.display = 'none';
           myplot4.style.display = '';
          document.getElementById('toggle_button').innerHTML = 'Show Distinct Peptide Per Million Observed Spectra Plot';
        } else {
           myplot3.style.display = '';
           myplot4.style.display = 'none';
          document.getElementById('toggle_button').innerHTML = 'Show Total Observed Spectra Plot';
        }
			}
		</script>
     <TABLE>
       <A NAME='<B>Sample Category</B>'></A><DIV CLASS="hoverabletitle"><B>Peptide Identification by Sample Category</B></DIV>
       <TR> 
        <TD><button type='button' id='toggle_button' onclick=toggle_plot()>Show Distinct Peptide Per Million Observed Spectra Plot</button>
           <div style="width: 80vw">
            &nbsp;<div id="plot_div3" style="display:none; width: 80vw"> </div>
            &nbsp;<div id="plot_div4"></div>
            <br><br>
           <div>
       </TD>
       </TR> 
    </TABLE>
		<script type="text/javascript" charset="utf-8">
      $plot_js
		</script>
  ~;
  return $chart;
  
}
sub displayProt_PTM_plotly{
   my $self = shift;
  my %args = @_;
  my $protein = $args{protein};
  my $data = $args{data};
  my $sequence = $args{seq},
  my $atlas_build_id = $args{atlas_build_id};

  my $dataTable = "var data=[";
  my $curpos   = "var curpos=[";
  my $totalObs = "var totalObs=[";
  my $maxnobs  = 0;
  my $seqlen = length($sequence);

  $sequence =~ s/\*.*//g;
  my @aas = split(//, $sequence);
  my $sep = '';
  foreach my $pos (0..$#aas){
    my $aa=$aas[$pos];
    my $nobs = '';
    if ($aa =~ /[STY]/){
      if (defined $data->{$protein}{$pos}){
        $dataTable .= "['$aa'";
        foreach my $c (qw (nP<.01 nP<.05 nP<.20 nP.2-.8 nP>.80 nP>.95 nP>.99 no-choice)){
          $dataTable .=",$data->{$protein}{$pos}{$c}";
          if ($data->{$protein}{$pos}{$c} > $maxnobs){
            $maxnobs = $data->{$protein}{$pos}{$c} ;
          }
          $nobs += $data->{$protein}{$pos}{$c};
        }
        $dataTable .= "],\n";
      } else {
        $dataTable .= "['$aa','','','','','','','',''],\n";
      }
    } else {
       $dataTable .= "['$aa','','','','','','','',''],\n";
   }
    $curpos .= "$sep";
    $curpos .=$pos+1 ;
    $totalObs .= "$sep";
    $totalObs .= $nobs;
    $sep = ",";
  }
  $dataTable =~ s/,$//;
  $dataTable =~ s/\n$//;
  $dataTable .= "];\n";
  $curpos   .= "];\n";
  $totalObs .= "];\n";

  my $trace = "var col1 =  getCol(data, 0);\n";
  my @colors = ('red','orange','purple','grey','#007eca','skyblue','green','black');
  my @names = ('<0.01','0.01-0.05','0.05-0.20','0.20-0.80','0.80-0.95','0.95-0.99','0.99-1.00','no-choice');
  foreach my $i(1..8){
    my $j=$i+1;
    my $k=$i-1;
    $trace .= qq~
    var col$j = getCol(data, $i);
    var trace$i = {
      x:curpos,
      y:col$j,
      marker: {color:'$colors[$k]'},
      type: 'bar',
      name:'$names[$k]'
    };
    ~;
  }

  my $plot_js = qq~
    function getCol(matrix, col){
       var column = [];
       for(var i=0; i<matrix.length; i++){
          column.push(matrix[i][col]);
       }
       return column;
    };
    $dataTable;
    $curpos;
    $totalObs;
		var myDiv = document.getElementById('protPTM_div');
    $trace
		var layout = {
			tickmode:'linear',
			tdick:50,
      annotations: [],
      hovermode:'x unified',
      barmode:'stack',
      yaxis:{title:'nObs',rangemode:'tozero'}
		};
   
		for ( var i = 0 ; i < $seqlen; i++ ) {
      if (totalObs[i] > 0){
				var result = {
					x: curpos[i], 
          y: totalObs[i], 
					text: '<b>'+totalObs[i]+'</b>',
					xanchor: 'center',
					yanchor: 'top',
          yshift:30,
					showarrow: false
				};
				layout.annotations.push(result);
     }
   }
	 plotdata= [trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8];
   Plotly.newPlot(myDiv, plotdata, layout, {scrollZoom: true});
	 var ticklabels=[];
	 for (var i in curpos){
			ticklabels[i]= curpos[i]+ '<br>'+col1[i];
	 }
	 myDiv.on('plotly_relayout',function(eventdata){  
		 //alert(JSON.stringify(eventdata));
			var xmax=eventdata['xaxis.range[1]'];
			var xmin=eventdata['xaxis.range[0]'];

			var n = xmin  < 0 ? 0:parseInt(xmin);
			var m = xmax  > $seqlen ? $seqlen-1 :parseInt(xmax);
			var ymax=0;
			var len = data.length;
			//get max obs
			var l =0;
			for (var i =n; i <=m; i++){
				var o = parseInt(totalObs[i]) || 0;
				l += o;
				if(o > ymax){	
					ymax=o;
				}    
			}
			ymax = ymax + 5;
			if ( xmax- xmin <= 40){
			 var update = {
				 xaxis:{
				 tickmode:'array',
					tickvals:curpos,
					ticktext:ticklabels,
					tickangle:0,
					range: [xmin, xmax], 
				 },
				 yaxis:{
					title:'nObs',
					range:[0,ymax]
				 }
			 };
				Plotly.relayout(myDiv, update);
			}
			if ( xmax- xmin > 40){	
				 if (xmin < 0){
					 xmin=0;
				 }
				 if (xmax> $seqlen){
					 xmax=$seqlen;
				 }
				 var update = {
					 xaxis:{
					 tickmode:'linear',
						 dtick:50,
						range: [xmin, xmax], 
					 },
					 yaxis:{
						title:'nObs',
						range:[0,ymax]
					 }
				 };
				 Plotly.relayout(myDiv, update);
			 }
			 if (eventdata['xaxis.autorange'] && eventdata['yaxis.autorange'] ){ 
				 var update = {
					 xaxis:{
					 tickmode:'linear',
					 dtick:50,
						yaxis:{title:'nObs',rangemode:'tozero'}
				 }
				};
				Plotly.relayout(myDiv, update);
			 }
	 });
  ~;

  my $chart = qq~
     <!-- Latest compiled and minified plotly.js JavaScript -->
     <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
     <div style="width: 80vw"><div id="protPTM_div"></div></div>
    <script type="text/javascript" charset="utf-8">
      $plot_js
    </script>
  ~;
  return $chart;

}
sub displayExperiment_contri_plotly{
  my $self = shift;
  my %args = @_;
  my $data_ref = $args{data_ref};
  my $tr = $args{tr};
  my $column_name_ref = $args{column_name_ref};
  my %cols =();
  my @sample_label = ();
  my (@cumpepx, @cumpepy,@idvpepy,@cumprotx, @cumproty,@idvproty);
  my $pre_cum_n_good_spectra;
  my $idx =0;
  foreach (@$column_name_ref){
    $cols{$_} = $idx;
    $idx++;
  } 
  $idx=0;
  foreach my $row(@$data_ref){
    my $n_good_spectra = $row->[$cols{"Spectra ID'd"}];
    my $n_distinct_peptides = $row->[$cols{"Distinct Peptides"}];
    my $cumulative_n_peptides = $row->[$cols{"Cumulative Peptides"}];
    my $n_canonical_proteins = $row->[$cols{"Distinct Canonical Proteins"}];
    my $cumulative_n_proteins = $row->[$cols{"Cumulative Canonical Proteins"}];
    my $sample_tag =  $row->[$cols{"Experiment Name"}];
    $sample_tag =~ s/.*sample_id=\d+\'>//;
    $sample_tag =~ s/<.*//;
    push @sample_label, ($sample_tag,'','');

		if ($idx ==0 ){
				push @cumpepx, 0;
			}else{
				push @cumpepx, $pre_cum_n_good_spectra;
        push @cumpepx, $pre_cum_n_good_spectra;
			}
			push @cumpepx, $pre_cum_n_good_spectra+$n_good_spectra;
      
			push @cumpepy, ($cumulative_n_peptides,$cumulative_n_peptides,0);
			push @idvpepy, ($n_distinct_peptides,$n_distinct_peptides,0);
			push @cumproty, ($cumulative_n_proteins,$cumulative_n_proteins,0);
			push @idvproty, ($n_canonical_proteins,$n_canonical_proteins,0);
			$pre_cum_n_good_spectra += $n_good_spectra;
			$idx++;
  }
  my $cumpepx_str = join(",", @cumpepx);
  my $cumpepy_str = join(",", @cumpepy);  
  my $idvpepy_str = join(",", @idvpepy);
  my $cumproty_str = join(",", @cumproty);
  my $idvproty_str = join(",", @idvproty);
  my $sample_label_str = join("','", @sample_label);
  my $total_spec = $cumpepx[$#cumpepx];
  my $width =  $total_spec < 5000000 ? 'width: 1000,' : '';
  my $plot_js = qq~
			l=['$sample_label_str']
			var cum = {
				x: [$cumpepx_str],
				y: [$cumpepy_str],
				fill: 'tonexty',
        fillcolor:'#1f77b4',
        marker: {
                  color: 'rgba(255,0,255,0)',
                }, 
				name:'cumulative_n_peptides',
				hovertext:l,
				haveron:"fills",
        hoverinfo:"x+y+text",
			};
			var idv = {
				x: [$cumpepx_str],
				y: [$idvpepy_str],
				fill: 'tonexty',
        fillcolor:'#ff7f0e',
				marker: {
					color: 'rgba(255,0,255,0)',
				},
				name:'n_distinct_peptides',
				hovertext:l,
        hoverinfo:"x+y+text",
				haveron:"fills"
			};

			var layout = {
				$width
				height: 800,
				font: {
					size: 18
				},
        legend:{x: 0.029,y: 1.1,font: { size: 12}},
				xaxis:{title:'Cumulative Number of MS/MS Spectra Identified'},
				yaxis:{title:'Number of Distinct Peptides'}
			};
			var data = [idv,cum];
			Plotly.newPlot('plot_div', data,layout);

      var cum2 = {
        x: [$cumpepx_str],
        y: [$cumproty_str],
        fill: 'tonexty',
        fillcolor:'#003f72',
        marker: {
                  color: 'rgba(255,0,255,0)',
                },
        name:'cumulative_n_proteins',
        hovertext:l,
        hoverinfo:"x+y+text",
        haveron:"fills"
      };
      var idv2 = {
        x: [$cumpepx_str],
        y: [$idvproty_str],
        fill: 'tonexty',
        fillcolor:'#b00',
        marker: {
                  color: 'rgba(255,0,255,0)',
                },
        name:'n_canonical_proteins',
        hovertext:l,
        hoverinfo:"x+y+text",
        haveron:"fills"
      };
      var layout = {
         $width
        height: 800,
        font: {
          size: 18
        },
        legend:{x: 0.029,y: 1.1,font: { size: 12}},
        xaxis:{title:'Cumulative Number of MS/MS Spectra Identified'},
        yaxis:{title:'Number of Distinct Proteins'}
      };

      var data2 = [idv2,cum2];
      Plotly.newPlot('plot_div2', data2,layout);
  ~;
  #print "$plot_js<BR>";

  my $chart = qq~
     <!-- Latest compiled and minified plotly.js JavaScript -->
     <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
     <div><p class=plot_caption><b>Plot below shows the number of peptides contributed by each experiment, and the cumulative number of distinct peptides for the build as of that experiment.</b></div> 
     <div style="width: 80vw"><div id="plot_div"></div><br></div></div>
     <div><p class=plot_caption><b>Plot below shows cumulative number of canonical proteins contributed by each experiment.</b><br>
     Height of blue bar is the number of proteins identified in experiment; 
     Height of red bar is the cumulatie number of proteins;<br> 
     Width of the bar (x-axis) shows the number of spectra identified (PSMs) above the threshold, for each experiment.</div>
     <div style="width: 80vw"><div id="plot_div2" ></div></div><br><br>
		<script type="text/javascript" charset="utf-8">
      $plot_js
		</script>
  ~;
  return $chart;
  
}
sub tableHeatMap{
  my $self = shift;
  my %args = @_;
  my $table_id = $args{table_id};
  my $total = $args{total} || 0;
  my $str = qq~
     <script type="text/javascript" src="$CGI_BASE_DIR/../usr/javascript/jquery/jquery.min.js"></script>
     <script type="text/javascript">
      jQuery.noConflict();
			jQuery(document).ready(function(){
				// Function to get the Max value in Array
					Array.sum = function( array ){
            var sum =0;
						for (let i in array){
              sum = sum+ array[i];
            }
            return sum;
					};

					// get all values
					var counts= jQuery('#$table_id  td:first-child+td').map(function() {
							if (jQuery(this).text() ){
								 return parseInt(jQuery(this).text());
							}else{
								 return 0;
							}
					}).get();
				// return max value
				var sum = Array.sum(counts);
        if ($total > 0){
           sum = $total;
        }
				// add classes to cells based on nearest 10 value
				jQuery('#$table_id  td:first-child + td').each(function(){
					var val = parseInt(jQuery(this).text());
					var pctval = parseInt((Math.round((val/sum)*100)).toFixed(0));
          var pctval2 = 100 - pctval ;
					clr = "linear-gradient(to right, #f46d69 " + pctval + "%, #EAEAEA " + pctval + "% " + pctval2 + "%)"; 
					jQuery(this).css("background-image", clr);
				});
			});
	 	</script>
  ~;
  return $str;
}

##################################################################################
sub plotly_barchart {
  my $self = shift;
  my %args = @_;
  my $all_data = $args{data};
  my $names = $args{names};
  my $divname = $args{divName} || '';
  my $xtitle = $args{xtitle} || '';
  my $ytitle = $args{ytitle} || '';
  my $title = $args{title} || '';
  my $dtick = $args{dtick} || 5;
  my $colors = $args{colors} ;
  
  my $plot_js = '';
  my $counter = 0;
  my @ts = ();
  foreach my $data (@$all_data){
		my @category = () ;
		my @cnt = () ;
		foreach my $row (@$data){
			push @category, $row->[0];
			push @cnt , $row->[1];
		}
		my $category_str = join(",", @category);
		my $cnt_str = join(",", @cnt);
    my $marker = '';
    if ($colors && $colors->[$counter]){
      $marker = "marker: {color: '$colors->[$counter]";
    }
		$plot_js .= qq~
				var t$counter = {
					x: [$category_str],
					y: [$cnt_str],
					name: '$names->[$counter]', 
					type: 'bar',	
				  $marker	
				};
			~;
     push @ts, "t$counter";
     $counter++;
   }
   my $ts_str = join(",", @ts);
   $plot_js .= qq~
				var data = [$ts_str];
				var layout = {
					barmode: 'group',
          font: {size: 18},
          legend:{x: 0.02,y: 1.25 ,font: { size: 12}},
          xaxis: {dtick:$dtick,title: '$xtitle',rangemode:'tozero'},
          yaxis: {title: '$ytitle'},
        };
        Plotly.newPlot('$divname', data, layout);
   ~;
  my $chart = qq~
    <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
		<script type="text/javascript" charset="utf-8">
      $plot_js
		</script>

  ~;
  return $chart;
}

sub get_dataset_url{
  my $self = shift;
  my $repository_ids = shift;
  my $url ='';
  my @ids = split(/[,;]/, $repository_ids);
  my %dataset_annotation = ();
	my $url = "http://proteomecentral.proteomexchange.org/api/autocomplete/v0.1/datasets";
	my $content  = get($url);
	$content =~ s/[\[\]"\s+\n\r]//g;
	foreach my $pxd (split(",", $content)){
		next if ($pxd eq '');
		$dataset_annotation{$pxd} =1;
	}
  $url='';

  foreach my $id(@ids){
    $id =~ s/\s+//g;
    if ($id =~ /PXD/){
      $url = "http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=";
    }elsif($id =~ /PASS/){
      $url = "http://www.peptideatlas.org/PASS/";
    }elsif($id =~ /^S\d+/){
      $url = "https://cptac-data-portal.georgetown.edu/study-summary/";
    }
    if ($url){
      $url= "<a href='$url$id' target='_blank'>$id</a>";
      if (defined $dataset_annotation{$id}){
        $url .= "<a href='http://proteomecentral.proteomexchange.org/devLM/annotation/view.html?dataset_id=$id' target='_blank'>&nbsp;[annot]</a>";
      }
      $url .= ",";
    }else{
      $url = "$id,";
    }
  }
  $url =~ s/,$//;
  return $url;
}


1;

__END__
###############################################################################
###############################################################################
###############################################################################

=head1 NAME

SBEAMS::WebInterface::HTMLPrinter - Perl extension for common HTML printing methods

=head1 SYNOPSIS

  Used as part of this system

    use SBEAMS::WebInterface;
    $adb = new SBEAMS::WebInterface;

    $adb->printPageHeader();

    $adb->printPageFooter();

    $adb->getGoBackButton();

    $adb->getGoBackButton();

=head1 DESCRIPTION

    This module is inherited by the SBEAMS::WebInterface module,
    although it can be used on its own.  Its main function 
    is to encapsulate common HTML printing routines used by
    this application.

=head1 METHODS

=item B<printPageHeader()>

    Prints the common HTML header used by all HTML pages generated 
    by theis application

=item B<printPageFooter()>

    Prints the common HTML footer used by all HTML pages generated 
    by this application

=item B<getGoBackButton()>

    Returns a form button, coded with javascript, so that when it 
    is clicked the user is returned to the previous page in the 
    browser history.

=head1 AUTHOR

Eric Deutsch <edeutsch@systemsbiology.org>

=head1 SEE ALSO

perl(1).

=cut


