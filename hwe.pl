#!"c:\Perl64\bin\perl.exe"	-w

use strict;
use warnings;
use Statistics::R;
 my $obs11=<STDIN>;
 my $obs12=<STDIN>;
 my $obs22=<STDIN>;
# Integrated Program for Allele Frequencies, Genotype Frequencies, HWE for diallelic markers
# rare homozygotes $obs_11
# common homozygotes $obs_22
my $obs_11;
my $obs_12;
my $obs_22;
if ($obs11 < $obs22)
  {
    $obs_11 = $obs11;
    $obs_22 = $obs22;
   }
else
  {
    $obs_11 = $obs22;
    $obs_22 = $obs11;
   }
 $obs_12 = $obs12;
 print $obs_11;
 print $obs_12;
 print $obs_22;
  
# Calculation for Allele Frequency

 #print " Allele Frequencies\n" . "<br/>";
# total number of genotypes
 my $nind = $obs_11 + $obs_12 + $obs_22;
 print $nind;
# print " nind: " . $nind . "<br/>";

#p value
#print " P value" . "<br/>";
my $allele_1=(($obs_11 * 2 + $obs_12)/(2 * $nind));
#print "allele_1:   " . $allele_1;

# q value
#print " Q value" . "<br/>";
 my $allele_2=(($obs_22 * 2 + $obs_12)/(2 * $nind));
#print " allele_2: " . $allele_2;

# Calculation of Genotype Frequency
#print " Genotype Frequency" . "<br/>";

my $genotype_1 = ($obs_11 / $nind);
my $genotype_2 = ($obs_12 / $nind);
my $genotype_3 = ($obs_22 / $nind);

print " Genotype_1: " . $genotype_1;
print " Genotype_2: " . $genotype_2;
print " Genotype_3: " . $genotype_3;

# Calculation of Expected Frequencies

#print " Expected Frequency" . "<br/>";
 my $expgenotype_1 = (($allele_1 * $allele_1) * $nind);
print " expgenotype_1:  " . $expgenotype_1;
#$deviation_1=$expgenotype_1 - $obs11;

 my $expgenotype_2 = ((2 * $allele_1 * $allele_2) * $nind);
print " expgenotype_2:  " . $expgenotype_2;
my $expgenotype_3 = (($allele_2 * $allele_2) * $nind);
print " expgenotype_3:  " . $expgenotype_3;

 my $enind=($expgenotype_1 + $expgenotype_2 + $expgenotype_3);
if(($expgenotype_1 < 5) || ($expgenotype_2 < 5) || ($expgenotype_3 < 5))
{ 
print "One or more of the expected genotypic frequencies is lower than 5";
print "The simulation approach should be used";
 }
else 
{
print "\n";
}
my $x1= ((($obs_11 - $expgenotype_1) * ($obs_11 - $expgenotype_1))/$expgenotype_1);
my $x2= ((($obs_12 - $expgenotype_2) * ($obs_12 - $expgenotype_2))/$expgenotype_2);
my $x3= ((($obs_22 - $expgenotype_3) * ($obs_22 - $expgenotype_3))/$expgenotype_3);
my $chisquare = ($x1 + $x2 + $x3); 

my $R = Statistics::R->new();
 $R->startR ;
 #$R->send('library(genetics);');
 $R->set('chisquare',$chisquare);
 $R->send('pp<-pchisq(chisquare,1,lower.tail=FALSE);');
 $R->send('print(pp);');
 $R->send('write.table(pp,"pp.txt");');
  #print "<pre>";
 my @pp = $R->get(('pp'));
 #print "P=",@pp,"\n";
 #print "</pre>";
  $R->stopR;
  
my $Amaybe = (($obs_12*$obs_12)/($obs_22*4));
my $pAfreq = ((($Amaybe + (0.5 * $obs_12)) / ($Amaybe + $obs_12 + $obs_22)*100))/100;
my $qAfreq = (100-($pAfreq*100))/100;
 my $Bmaybe = sqrt($obs_11*$obs_22*4);
my $pBfreq = ((($obs_11 + (0.5 * $Bmaybe)) / ($obs_11 + $Bmaybe + $obs_22)*100))/100;
 my $qBfreq = (100-($pBfreq*100))/100;
 my $Cmaybe = (($obs_12*$obs_12)/($obs_11*4));
my $pCfreq = ((($obs_11 + (0.5 * $obs_12)) / ($obs_11 + $obs_12 + $Cmaybe)*100))/100;
my $qCfreq = (100-($pCfreq*100))/100;
my $nosamples = ($obs_11 + $obs_12 + $obs_22);
#my $sumofchi = 0;
#$sumofchi = $x1 + $x2 + $x3;
#my $significance = "";
#$significance = "<br>for likelihoods of calculated <i>X</i><sup>2</sup> value see below. ";
print "Chisquare=###.#### with 1 degree of freedom";

if ($chisquare < 3.841) 
{print "Chisquare is not significant at a 0.05 level";
}
if (($chisquare >= 3.8410) && ($chisquare < 6.634))
{ print "0.05 > p > 0.01";}
if(($chisquare >= 6.634) && ($chisquare < 10.827))
{print "0.01 > p > 0.001";
}
if ($chisquare >= 10.827) 
{print "p < 0.001";}
print "\n";

if ($chisquare > 3.84)

{
 my $chisquare1="It is significant (as the value is >3.84)";
#print "\nIt is significant" . "<br/>";
}
else
{
 my $chisquare1="It is not significant";
#print "\n It is not significant" . "<br/>";
}
	print "Common Homozygotes\n";
	print "Expected" . $x1. "\n";  
 	print "observed" . $obs11 . "\n";
	print "Heterozygotes\n";
	print "Expected" . $x2 . "\n"; 
	print "Observed" . $obs12 . "\n";
	print "Rare Homozygotes\n";
	print "Expected" . $x3 . "\n";
	print "Observed" . $obs22 . "\n";
	print "\np allele freq = ";
	print "$allele_1" . "\n";
	print "\nq allele freq = ";
	print "$allele_2" . "\n";
	print "\nchi-square = ";
	print "$chisquare" . "\n";
	print "\nP value for chisquare=";
	print "@pp";
	print "\n\n";
	print "$Amaybe";
	print "$obs12";
	print "$obs22";
	print "$pAfreq";
	print "$qAfreq";
	print "$obs11";
	print "$Bmaybe";
	print "$obs22";
	print "$pBfreq";
	print "$qBfreq";
	print "$obs11";
	print "$obs12";
	print "$Cmaybe";
	print "$pCfreq";
	print "$qCfreq";	
	