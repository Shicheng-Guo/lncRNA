use Cwd;
chdir getcwd;
my @file=glob("NA-*");
foreach my $file(@file){
my (undef,$newname)=split /\-/,$file;
system("cp $file $newname");
}

