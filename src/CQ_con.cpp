PQ_mat_struct CQ_con_cal(PQ_mat_struct &P,PQ_mat_struct &Q,int &coe)
{
	PQ_mat_struct temp;
	struct node<int *,int> *PQ_root=NULL;
	
	temp.M=P.M-Q.M;
	temp.pi=(P.pi+Q.pi)%2;
	if(Q.J==0)
	{
		temp.J=P.J;
	}
	else
	{
		temp.J=-1;
	}
	temp.PQ_Jpi_index=-1;
	for(int i=1;i<=P.val[0];i++)
	{
		for(int j=1;j<=Q.val[0];j++)
		{
			int size_P=P.ser[i][0];
			int size_Q=Q.ser[j][0];
			int size=P.ser[i][0]+Q.ser[j][0];
			int ser[size+1];
			ser[0]=size;
			int val=P.val[i]*Q.val[j];
			memcpy(ser+1,P.ser[i]+1,sizeof(int)*size_P);
			if(P.ser[i][size_P-1]==Q.ser[j][1])
			{
				ser[0]-=2;
				size-=2;
				ser[size_P]+=Q.ser[j][2];
				memcpy(ser+size_P+1,Q.ser[j]+3,sizeof(int)*(size_Q-2));
			}
			else
			{
				memcpy(ser+size_P+1,Q.ser[j]+1,sizeof(int)*size_Q);
			}
			PQ_tree_search_add(PQ_root,ser,val,PQ_cmp);
		}
	}
	int num_ser=P.val[0]*Q.val[0];
	temp.val=new int [num_ser+1];
	temp.ser=new int *[num_ser+1];
	temp.ser[0]=NULL;
	int index_ser=1;
	PQ_tree_op(PQ_root,temp.ser,temp.val,index_ser);
	temp.val[0]=index_ser-1;
	if(temp.val[0]!=0)
	{
		coe=abs(temp.val[1]);
		for(int i=2;i<=temp.val[0];i++)
		{
			coe=gcb(coe,temp.val[i]);
		}
		if(temp.val[1]<0)
		{
			coe=-coe;
		}
		for(int i=1;i<=temp.val[0];i++)
		{
			temp.val[i]/=coe;
		}
	}
	else
	{
		coe=0;
	}
	return temp;
}