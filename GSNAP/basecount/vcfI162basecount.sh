

cat $1 | awk '{print $1, $2, $4, $5, $8;}' | perl -e 'while(<STDIN>){chomp; my ($chr, $pos, $ref, $alt, $info) = split; if($info =~ m/I16=(\d+,\d+,\d+,\d+)/){ $dp4=$1; @dp4arr=split(",", $dp4); $nref = $dp4arr[0] + $dp4arr[1]; $nalt = $dp4arr[2] + $dp4arr[3]; $tot = $nref + $nalt; if($nalt != 0){$frac = $nalt / $tot}else{$frac = 0;}; print join("\t", ($chr, $pos, $ref, $alt, $tot, $nref, $nalt, $frac)) . "\n";}}'
