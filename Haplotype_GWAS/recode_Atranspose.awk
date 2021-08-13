{
	for(i = 7; i <=  NF; i++)
	gsub(0,-1,$i);
	gsub(1,0,$i);
	gsub(2,1,$i);	
}
END {print $0}


