double Rout(int i,int j,int lambda,vector<nlj_struct> orbit)//lambda does not need *2
{
	int orbit_num=orbit.size();
	if((orbit[i].l+orbit[j].l+lambda)%2!=0)
	{
		return 0;
	}
	int nr[orbit_num];
	for(int k=0;k<=orbit_num-1;k++)
	{
		nr[k]=(orbit[k].N-orbit[k].l)/2;
	}
	double res=sgn(nr[j]-nr[i])*sqrt(fac_out(nr[j])/fac2_out(2*nr[j]+2*orbit[j].l+1)*fac_out(nr[i])/fac2_out(2*nr[i]+2*orbit[i].l+1)*pow(2,nr[j]+nr[i]-lambda));
	res=res*fac_out((orbit[i].l-orbit[j].l+lambda)/2)*fac_out((orbit[j].l-orbit[i].l+lambda)/2);
	double sum_q=0;
	double delta_q;
	int q,q_min,q_max;
	q_min=0;
	if(nr[i]-(orbit[j].l-orbit[i].l+lambda)/2>0)
	{
		q_min=nr[i]-(orbit[j].l-orbit[i].l+lambda)/2;
	}
	if(nr[j]-(orbit[i].l-orbit[j].l+lambda)/2>q_min)
	{
		q_min=nr[j]-(orbit[i].l-orbit[j].l+lambda)/2;
	}
	q_max=nr[j];
	if(nr[i]<q_max)
	{
		q_max=nr[i];
	}
	for(q=q_min;q<=q_max;q++)
	{
		delta_q=fac2_out(orbit[j].l+orbit[i].l+lambda+2*q+1)/pow(2,q)/fac_out(q)/fac_out(nr[j]-q)/fac_out(nr[i]-q)/fac_out(q+(orbit[j].l-orbit[i].l+lambda)/2-nr[i])/fac_out(q+(orbit[i].l-orbit[j].l+lambda)/2-nr[j]);
		sum_q=sum_q+delta_q;
	}
	res=res*sum_q;
	return(res);
}

double Y_ring(int i,int j,int lambda,vector<nlj_struct> orbit)//lambda does not need *2
{
	int orbit_num=orbit.size();
	if((orbit[i].l+orbit[j].l+lambda)%2!=0)
	{
		return(0);
	}
	else
	{
		double res=sgn(orbit[i].l+orbit[j].l+(orbit[j].j-orbit[i].j)/2)*sqrt((2*lambda+1)*(orbit[j].j+1)*(orbit[i].j+1)/4.0/PI);
		int jser[6];
		jser[0]=orbit[i].j;
		jser[1]=lambda*2;
		jser[2]=orbit[j].j;
		jser[3]=-1;
		jser[4]=0;
		jser[5]=1;
		res=res*sgn((orbit[i].j-1)/2)*threej_out(jser);
		res=res/sqrt(orbit[i].j+1);
		return(res);
	}
}
int J_pm_index[2];

void J_pm_out(int norp)
{
	double *q=new double [sp_nljm_dim2[norp]];
	memset(q,0,sizeof(double)*sp_nljm_dim2[norp]);
	for(int i=0;i<sp_nljm[norp].size();i++)
	{
		for(int j=0;j<sp_nljm[norp].size();j++)
		{
			if(sp_nljm[norp][i].m==sp_nljm[norp][j].m+2&&sp_nljm[norp][i].nlj==sp_nljm[norp][j].nlj)
			{
				q[i+j*sp_nljm[norp].size()]=sqrt((sp_nljm[norp][i].j-sp_nljm[norp][j].m)*(sp_nljm[norp][i].j+sp_nljm[norp][j].m+2))/2.0;
			}
		}
	}
    PQ_M_struct temp;
	temp.J=2;
	temp.M=2;
	temp.pi=0;
	temp.PQ_Jpi_index=-1;
	
	double *qT=new double [sp_nljm_dim2[norp]];
	memcpy(qT,q,sizeof(double)*sp_nljm_dim2[norp]);
	mkl_dimatcopy('c','t',sp_nljm[norp].size(),sp_nljm[norp].size(),1,qT,sp_nljm[norp].size(),sp_nljm[norp].size());
	temp.yq=q;
	temp.yq_mM=qT;
	Q_M[norp].emplace_back(temp);
	J_pm_index[norp]=Q_M[norp].size()-1+INT_MAX/2;
}