int q1 = 1000/100 *20;

if( !(i/double(q1)-i/q1) )
{
	ofstream * file = new ofstream("out.m",ios::app);

	(*file)<<"xs=[xs "<<x[16580-1]<<"];"<<endl ;
	
	file->close();
	delete file;
}

int q = 100000 / 20 ; // 100;

if( !(i/double(q)-i/q) )
{
	char str1[80]="outpt/out";
	char str3[80]=".m";

	strncat(str1,itoa(i/q,10),10);
	strncat(str1,str3,10);

	ofstream * file = new ofstream(str1);

	// file->precision(20);

	(*file)<<" ";
	for(j=0; j<NumPoints; j++)
	{
		(*file)<<x[j]<<" ";
	}
	(*file)<<" "<<endl;

	(*file)<<" ";
	for(j=0; j<NumPoints; j++)
	{
		(*file)<<y[j]<<" ";
	}
	(*file)<<" "<<endl;

	file->close();
	delete file;
}
