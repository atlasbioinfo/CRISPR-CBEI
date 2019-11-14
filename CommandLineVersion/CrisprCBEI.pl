use strict;use warnings;
our %baseT = (
    "A"=> "T",
    "T"=> "A",
    "C"=> "G",
    "G"=> "C",
    "U"=> "A",
    "a"=> "t",
    "t"=> "a",
    "c"=> "g",
    "g"=> "c",
    "u"=> "a",
    "N"=> "N",
    "n"=> "n",
    "R"=> "Y",
    "Y"=> "R",
    "M"=> "K",
    "K"=> "M",
    "W"=> "W",
    "S"=> "S",
    "H"=> "D",
    "D"=> "H",
    "B"=> "V",
    "V"=> "B",
    "r"=> "y",
    "y"=> "r",
    "m"=> "k",
    "k"=> "m",
    "w"=> "w",
    "s"=> "s",
    "h"=> "d",
    "d"=> "h",
    "b"=> "v",
    "v"=> "b"
);

open F,"$ARGV[0]";
my $count=0;
my $pam="NGG";
$pam=~s/\s|\r|\t|\n//g;
my $pamr=&tranStrand($pam);
my $reg=&getPAMreg($pam);
my $regr=&getPAMreg($pamr);
my $spacer=20;
my $ebeg=4;
my $eend=8;
my $dire=3;
{
	local $/=">";
	while (<F>){
		chomp;
		next if length($_)<2;
		my @sp=split(/\n/);
		my @title=split(/\s/,$sp[0]);
		my $seq="";
		for my $i(1..$#sp){
			$sp[$i]=~s/\s|\r|\n|\t//g;
			$seq.=$sp[$i]
		}
		my @codon=$seq=~/(\w{3})/g;
		for my $j(0..$#codon){
			if ($codon[$j]=~/(CAA|CAG|CGA)/i){
				my $tpos=$j*3;
				if ($dire==5){
					my $tb=$tpos+$spacer-$eend+1;
					my $tl=$eend-$ebeg+length($pam);
					next if ($tb<0 || $tb+$tl-1>length($seq));
					my $tarseq=substr($seq,$tb,$tl);
					while ($tarseq=~m/(?<=$reg)/ig){
						my $ttpos=pos($tarseq);
						next if $tb+$ttpos-1>length($seq);
						print "$title[0]\tPlus\t";#Fasta title
						print substr($seq,$tb+$ttpos-length($pam)-$spacer,$spacer),"\t";
						print $tb+$ttpos-length($pam)-$spacer+1,"-",$tb+$ttpos-length($pam),"\t";
						print substr($seq,$tb+$ttpos-length($pam),length($pam)),"\t"; #PAM seq
						print $tb+1+$ttpos-length($pam),"-",$tb+$ttpos,"\t"; #PAM pos
						print "\n";
					}
				}else{
					my $tb=$tpos-$eend-length($pam)+1;
					my $tl=$eend-$ebeg+length($pam);
					next if ($tb<0 || $tb+$tl-1>length($seq));
					my $tarseq=substr($seq,$tb,$tl);
#					print "$tpos\t$tb\t$tl\t",$tarseq,"\n";
					while ($tarseq=~m/(?<=$reg)/ig){
						my $ttpos=pos($tarseq);
						next if $tb+$ttpos+$spacer-1>length($seq);
						print "$title[0]\tPlus\t";#Fasta title
						print &tranStrand(substr($seq,$tb+$ttpos,$spacer)),"\t";
						print $tb+$ttpos+1,"-",$tb+$ttpos+$spacer,"\t";
						print substr($seq,$tb+$ttpos-length($pam),length($pam)),"\t"; #PAM seq
						print $tb+$ttpos-length($pam)+1,"-",$tb+$ttpos,"\n"; #PAM pos
					}
					
				}
			}elsif($codon[$j]=~/(TGG)/i){
				my $tpos=$j*3;
				if ($dire==5){
					my $tb=$tpos+$ebeg-$spacer-length($pam)+1;
					my $tl=$eend-$ebeg+1+length($pam);
					next if ($tb<0 || $tb+$tl-1>length($seq));
					my $tarseq=substr($seq,$tb,$tl);
					while ($tarseq=~m/(?<=$regr)/ig){
						my $ttpos=pos($tarseq);
						next if $tb+$ttpos+$spacer-1>length($seq);
						print "$title[0]\tMinus\t";#Fasta title
						print &tranStrand(substr($seq,$tb+$ttpos,$spacer)),"\t";
						print $tb+$ttpos+1,"-",$tb+$ttpos+$spacer,"\t";
						print &tranStrand(substr($seq,$tb+$ttpos-length($pam),length($pam))),"\t"; #PAM seq
						print $tb+$ttpos-length($pam)+1,"-",$tb+$ttpos,"\n"; #PAM pos
					}
				}else{
					my $tb=$tpos+$ebeg+1;
					my $tl=$eend-$ebeg+length($pam)+1;
					next if ($tb<0 || $tb+$tl-1>length($seq));
					my $tarseq=substr($seq,$tb,$tl);
					while ($tarseq=~m/(?<=$regr)/ig){
						my $ttpos=pos($tarseq);
						next if $tb+$ttpos-1>length($seq);
						print "$title[0]\tMinus\t";#Fasta title
						print substr($seq,$tb+$ttpos-length($pam)-$spacer,$spacer),"\t";
						print $tb+$ttpos-length($pam)-$spacer+1,"-",$tb+$ttpos-length($pam),"\t";
						print &tranStrand(substr($seq,$tb+$ttpos-length($pam),length($pam))),"\t"; #PAM seq
						print $tb+1+$ttpos-length($pam),"-",$tb+$ttpos,"\t"; #PAM pos
						print "\n";
					}
					
				}
			}
		}
	}
}


sub getPAMreg(){
	my $pam=shift(@_);
	my @in=split(//,$pam);
	my $reg="";
	for my $i(0..$#in){
	    if ($in[$i]=~/N/g){
	        $reg.="\\w";
	    }elsif($in[$i]=~/[A,T,C,G]/g){
	        $reg.=$in[$i];
	    }elsif($in[$i]=~/R/g){
	        $reg.="[A,G,R]";
	    }elsif($in[$i]=~/Y/g){
	        $reg.="[C,T,Y]";
	    }elsif($in[$i]=~/M/g){
	        $reg.="[A,C,M]";
	    }elsif($in[$i]=~/K/g){
	        $reg.="[G,T,K]";
	    }elsif($in[$i]=~/S/g){
	        $reg.="[G,C,S]";
	    }elsif($in[$i]=~/W/g){
        	$reg.="[A,T,W]";
	    }elsif($in[$i]=~/H/g){
        	$reg.="[A,T,C,W,M,Y,H]";
	    }elsif($in[$i]=~/B/g){
        	$reg.="[G,T,C,K,S,Y,B]";
	    }elsif($in[$i]=~/V/g){
	        $reg.="[G,A,C,R,S,M,V]";
	    }elsif($in[$i]=~/D/g){
        	$reg.="[G,A,T,R,K,W,D]";
	    }else{
        	$reg.=$in[$i];
	    }
	}
	return $reg;	
}

sub tranStrand(){
	my $seq=shift(@_);
	my @sp=split(//,$seq);
	@sp=reverse(@sp);
	my $new="";
	for my $i(0..$#sp){
		$new.=$baseT{$sp[$i]};
	}
	return $new;
}
