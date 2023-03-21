#!/usr/bin/perl

# Author: Kevin Wiehe, Duke Human Vaccine Institute, 2016

# Parse MPEx results *_MPExResults.txt file into two text files for each sequence:
# - <sequence_name>.residue.txt (every line is a sequence position and the amino acid residue at this position)
# - <sequence_name>.dwif.txt (every line is a sequence position and the dwif score for the corresponding amino acid)

if (@ARGV==0){print "$0 [fasta] [MPEx results file]\n"; exit(1);}

open(FILE, $ARGV[0]);
@lines=<FILE>;
close(FILE);

for($i=0; $i<@lines; $i+=2)
{
    chomp $lines[$i];
    chomp $lines[$i+1];
    if ($lines[$i]!~/^>/){print STDERR "1 seq per header error\n"; exit(1);}
    ($name)=$lines[$i]=~/^>(.*?)$/;
    push(@names, $name);
    push(@sequences, $lines[$i+1]);
}

open(FILE, $ARGV[1]);
@lines=<FILE>;
close(FILE);
$NAME_COUNTER=-1;
for($i=0; $i<@lines; $i++)
{
    chomp $lines[$i];
    @a=split(/\t/, $lines[$i]);
#    print "$lines[$i]\n"; 
    if ($a[0]=~/^\"Position/){$NAME_COUNTER++; }
    else
    {
	$hash{$names[$NAME_COUNTER]}{$a[0]*1}{"res"}=$a[1];
	$hash{$names[$NAME_COUNTER]}{$a[0]*1}{"dG"}=$a[3];
	$hash{$names[$NAME_COUNTER]}{$a[0]*1}{"win"}=$a[4];
	$hash{$names[$NAME_COUNTER]}{$a[0]*1}{"swin"}=$a[5];
#	print "$NAME_COUNTER\t$names[$NAME_COUNTER]\n"; 
    }
}

foreach $name (@names)
{
    open(OUT, ">$name.smooth_hydroph.txt");
    for($pos=1; $pos<=200; $pos++)
    {
	if ($hash{$name}{$pos}{swin}=~/\d/){
	    print OUT "$pos\t$hash{$name}{$pos}{swin}\n"; }
    }
    close(OUT);

    open(OUT, ">$name.residue.txt");
    for($pos=1; $pos<=200; $pos++)
    {
	if (defined($hash{$name}{$pos}{res}))
	{ print OUT "$pos\t$hash{$name}{$pos}{res}\n"; }
    }
    close(OUT);

    open(OUT, ">$name.dwif.txt");
    for($pos=1; $pos<=200; $pos++)
    {
	if ($hash{$name}{$pos}{dG}=~/\d/)
	{ print OUT "$pos\t$hash{$name}{$pos}{dG}\n"; }
    }
    close(OUT);
    
}
