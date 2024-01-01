#!perl

use strict;
use warnings;
use Getopt::Long;
use MaterialsScript qw(:all);

use IO::Handle;

use Data::Dumper qw(Dumper);
use Time::HiRes qw(gettimeofday);

# To make directroy.
use File::Path qw(make_path);

# Set simulation pahts.
my $importDir = "E:/1_RODMOF/cand1-gen/1_beforeMS";
my $exportDir = $importDir ."_optimized";
my $energyOut = $exportDir . "/energy.txt";

make_path $exportDir;

my @cifPaths = glob($importDir . "/*.cif");
#my $runListFile = $Documents{"runlist.txt"};
my $dirName = "test";

my $start = gettimeofday();

#foreach my $docName (@{$runListFile -> Lines})
foreach my $cifPath (@cifPaths)
{
    # Remove parent folders.
    my @splitted = split/\//,$cifPath;
    # Get cif name.
    my $docName = $splitted[-1];
    
    # SKIP existing file.
    if (-e "$exportDir/$docName") {
    	print "SKIP existing file. ", "$exportDir/$docName \n";
        next;
    }
    
    # Remove suffix.
    my $rev = reverse $docName;
    my $txt = substr($rev, 4);
    $docName = reverse $txt;
    
    print "Optimization of $docName starts.\n";
    
    # Load MS document.
    #$docName =~ s/[\n\r]//g;
    my $doc = "";
    my $cif = "$importDir/$docName.cif";
    eval {
    	$doc = Documents -> Import($cif);
    	1;
    } or do {
    	print "Import error. $cif, SKIP. \n";
    	next;
    };
    $doc -> Close;

    # Load cif data (bond type, pbc info)
    my %bond_types;
    my %bond_pbc;
    open(my $fh, "<", $cif) || die "Couldn't open '".$cif."' because: ".$!;
    while(my $line = <$fh>){
       $line =~ s/\s+$//;
       last if ($line eq "_ccdc_geom_bond_type");
    }

    while(my $line = <$fh>){
       #print $line;
        my @sline = split /\s+/, $line;
        #print Dumper \@sline;
        $bond_types{$sline[0]}{$sline[1]} = $sline[4];
        $bond_pbc{$sline[0]}{$sline[1]} = $sline[3];
        $bond_types{$sline[1]}{$sline[0]} = $sline[4];
        $bond_pbc{$sline[1]}{$sline[0]} = $sline[3];
    }
    close $fh;
    
    #print Dumper \%bond_types;
    
    # Matching bonds (cif & MS document)
    my $mserror = "";
    my $msbonds = $doc->UnitCell->Bonds;
    foreach my $bond (@$msbonds) {
       my $atom1 = $bond->Atom1->Name;
       my $atom2 = $bond->Atom2->Name;
       my $mschemtype = substr($bond->ChemicalType, 0, 1);
       my $cifchemtype = $bond_types{$atom1}{$atom2};
       
       if (not $cifchemtype) {
           print "Empty string in cifchemtype, Neglect the structure.\n";
           $mserror = 1;
       }
       
       if ($mserror) {
           last;
       }
       
       if ($mschemtype ne $cifchemtype) {
           my $newbondtype = "";
           if ($cifchemtype eq "A"){
              # Not sure about this.
              $newbondtype = "Partial double";
           }
           elsif ($cifchemtype eq "S"){
              $newbondtype = "Single";
           }
           elsif ($cifchemtype eq "D"){
              $newbondtype = "Double";
           }
           elsif ($cifchemtype eq "T"){
              $newbondtype = "Triple";
           }
           else{
              printf "$docName Bond Type Error, ($cifchemtype). Neglect the bond.\n";
              printf "%s-%s : %s %s\n", $atom1, $atom2, $mschemtype, $cifchemtype;
              next;
           }
           $bond->BondType=$newbondtype;
           #printf "%s, %s\n\n", $bond->BondType, $bond->ChemicalType;
       }
    }	
    
    if ($mserror) {
        $doc -> Close;
        $doc -> Delete;
        next;
    }
    
    my $atoms = $doc->UnitCell->Atoms;
    foreach my $atom (@$atoms) {
        my $atomType = $atom->ElementSymbol;
	if ($atomType eq 'He')
	    {
	    $atom->Delete;
	    }
        }

    my $lenB = "LengthB";	# Axis which is compressed
    my $lenC = "LengthC";	# Axis which is stretched
    my $lenA = "LengthA";	# Axis which is stretched

    my $newA = $doc->SymmetryDefinition->LengthA;
    my $newB = $doc->SymmetryDefinition->LengthB;
    my $newC = $doc->SymmetryDefinition->LengthC;
    $newA = $newA * 1.2;
    $newB = $newB * 1.2;
    $newC = $newC * 1.2;
    
    $doc -> UnitCell -> Lattice3D -> $lenA = $newA;
    $doc -> UnitCell -> Lattice3D -> $lenB = $newB;
    $doc -> UnitCell -> Lattice3D -> $lenC = $newC;
    
    eval {
        # Optimize and Write energy outputs
        my $results = Modules->Forcite->GeometryOptimization->Run($doc, Settings(
   	    	          OptimizationAlgorithm => "Smart", 
                          Quality => "Medium",
                          ChargeAssignment => 'Use current',
                          AssignBondOrder => 'No',     
                          OptimizeCell => 'No',
                          WriteLevel => "Silent"));

        my $results = Modules->Forcite->GeometryOptimization->Run($doc, Settings(
   	    	          OptimizationAlgorithm => "Smart", 
                          Quality => "Fine",
                          ChargeAssignment => 'Use current',
                          AssignBondOrder => 'No',     
                          OptimizeCell => 'Yes',
                          WriteLevel => "Silent"));

        open(my $fhout, ">>", $energyOut) || die "Couldn't open '".$energyOut."' because: ".$!;
        printf $fhout "$docName %12.6f\n", $doc->PotentialEnergy;
        close $fhout;
        $doc -> Export("$exportDir/$docName.cif");
        1;
    } or do {
        print "Optimization fails, SKIP. \n"; 
    };
    $doc -> Close;
    $doc -> Delete;
}
my $end = gettimeofday();
my $runtime = $end - $start;

print $runtime, " sec";

