PQ_mat_struct CC_con_cal(PQ_mat_struct &P1,PQ_mat_struct &P2)
{
	PQ_mat_struct temp;
	temp.M=P2.M-P1.M;
	temp.pi=(P1.pi+P2.pi)%2;
	temp.PQ_Jpi_index=-1;
	if(P1.J==0||P2.J==0)
	{
		temp.J=P1.J+P2.J;
	}
	else
	{
		temp.J=-1;
	}
	int num=P2.val[0]*P1.val[0];
	temp.val=new int [num+1];
	temp.ser=new int *[num+1];
	temp.val[0]=num;
	temp.ser[0]=NULL;
	num=1;
	for(int i=1;i<=P2.val[0];i++)
	{
		for(int j=1;j<=P1.val[0];j++)
		{
			int size_P2=P2.ser[i][0];
			int size_P1=P1.ser[j][0];
			int size=size_P2+size_P1;
			int ser[size+1];
			int val=-P2.val[i]*P1.val[j];
			ser[0]=size;
			memcpy(ser+1,P2.ser[i]+1,sizeof(int)*size_P2);
			memcpy(ser+size_P2+1,P1.ser[j]+1,sizeof(int)*size_P1);
			temp.ser[num]=new int [size+1];
			memcpy(temp.ser[num],ser,sizeof(int)*(size+1));
			temp.val[num]=val;
			num++;
		}
	}
	return temp;
}