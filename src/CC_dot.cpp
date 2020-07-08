double C_dot_sum(node<int*,int> *&dot,int &norp)
{
	if(dot!=NULL)
	{
        double res=C_dot_sum(dot->l,norp);
        if(dot->v!=0)
        {
            res+=dot->v*C_tr_cal(dot->k,norp);
        }
        delete [] dot->k;
		res+=C_dot_sum(dot->r,norp);
        free(dot);
        return res;
	}
    else
    {
        return 0;
    }
    
}
double CC_dot_cal(struct PQ_mat_struct & P1,struct PQ_mat_struct & P2,int norp)
{
    if(P1.M!=P2.M)
    {
        return 0;
    }
    if(P1.pi!=P2.pi)
    {
        return 0;
    }
    if(P1.J!=-1&&P2.J!=-1&&P1.J!=P2.J)
    {
        return 0;
    }
    //stop here***********************
    struct node<int*,int> * dot=NULL;
	for(int i=1;i<=P1.val[0];i++)
	{
		for(int j=1;j<=P2.val[0];j++)
		{
            int size1=P1.ser[i][0];
            int size2=P2.ser[j][0];
            int size=size1+size2;
			int val=P1.val[i]*P2.val[j];
            int ser[size+1];
            ser[0]=size;
            memcpy(ser+1,P1.ser[i]+1,sizeof(int)*(P1.ser[i][0]));
			if(P1.ser[i][P1.ser[i][0]-1]==P2.ser[j][1]&&(abs(P2.ser[j][1])<INT_MAX/4||abs(P2.ser[j][1])>INT_MAX/2))
			{
                ser[size1]+=P2.ser[j][2];
                memcpy(ser+size1+1,P2.ser[j]+3,sizeof(int)*(size2-2));
                ser[0]-=2;
                size-=2;
			}
			else
			{
                memcpy(ser+size1+1,P2.ser[j]+1,sizeof(int)*size2);
			}
            if(ser[1]>ser[size-1])
            {
                int valT;
                int serT[size+1];
                qT(val,ser,valT,serT);
                dot_tree_search_add(dot,serT,valT,dot_cmp);
            }
            else
            {
                dot_tree_search_add(dot,ser,val,dot_cmp);
            }
		}
	}
    return C_dot_sum(dot,norp);
}