struct PQ_Jpi_struct
{
	int J;
	int pi;
	double *yq;
};
vector<PQ_Jpi_struct> P_Jpi[2];
vector<PQ_Jpi_struct> Q_Jpi[2];

struct PQ_M_struct
{
    int J;
    int pi;
    int PQ_Jpi_index;
	int M;
	double *yq;
    double *yq_mM;
};

vector<PQ_M_struct> P_M[2];
vector<PQ_M_struct> Q_M[2];

struct PQ_mat_struct
{
    int J;
    int pi;
	int PQ_Jpi_index;
	int M;
	int *val;
	int **ser;
};
inline void PQ0_add_init()
//this is set to let the Q0^T=Q0 case always stands.
{
	for(int norp=0;norp<2;norp++)
	{
		PQ_M_struct Q0;
		Q0.J=0;
		Q0.M=0;
		Q0.pi=0;
		Q0.PQ_Jpi_index=-1;
		Q0.yq=NULL;
		Q0.yq_mM=NULL;
		Q_M[norp].emplace_back(Q0);
		P_M[norp].emplace_back(Q0);
	}
}

inline void P_anti_symmetric(double *y,int norp)
{
	double P_temp[sp_nljm_dim2[norp]];
	memcpy(P_temp,y,sizeof(double)*sp_nljm_dim2[norp]);
	mkl_dimatcopy('c','t',sp_nljm[norp].size(),sp_nljm_dim2[norp],-1,y,sp_nljm[norp].size(),sp_nljm[norp].size());
	double alpha=1;
	int inc=1;
	daxpy(sp_nljm_dim2+norp,&alpha,P_temp,&inc,y,&inc);
}

int PJpi_to_PM(int P_i,int M,int norp)//the order of yq in PJpi follows the C stander, while P_M is follows Fortran stander.
{
	double *yq=new double [sp_nljm_dim2[norp]];
	memset(yq,0,sizeof(double)*sp_nljm_dim2[norp]);
	for(vector<nljm_struct>::iterator i=sp_nljm[norp].begin();i!=sp_nljm[norp].end();i++)
	{
		for(vector<nljm_struct>::iterator j=sp_nljm[norp].begin();j!=sp_nljm[norp].end();j++)
		{
			if(i->m+j->m==M&&P_Jpi[norp][P_i].J>=abs(i->j-j->j)&&P_Jpi[norp][P_i].J<=i->j+j->j&&(i->l+j->l)%2==P_Jpi[norp][P_i].pi)
			{
				yq[distance(sp_nljm[norp].begin(),i)+distance(sp_nljm[norp].begin(),j)*sp_nljm[norp].size()]+=P_Jpi[norp][P_i].yq[i->nlj*sp_nlj[norp].size()+j->nlj]*cg_cal(i->j,i->m,j->j,j->m,P_Jpi[norp][P_i].J,M);
				// if(P_Jpi[norp][P_i].J==0)
				// {
				// 	cout<<i->j<<' '<<i->m<<' '<<j->j<<' '<<j->m<<' '<<P_Jpi[norp][P_i].J<<' '<<M<<' '<<cg_cal(i->j,i->m,j->j,j->m,P_Jpi[norp][P_i].J,M)<<' '<<P_Jpi[norp][P_i].yq[i->nlj*sp_nlj[norp].size()+j->nlj]<<endl;
				// }
			}
		}
	}
	//cin>>wat;
	P_anti_symmetric(yq,norp);
	// if(P_Jpi[norp][P_i].J==0)
	// {
	// 	for(vector<nljm_struct>::iterator i=sp_nljm[norp].begin();i!=sp_nljm[norp].end();i++)
	// 	{
	// 		for(vector<nljm_struct>::iterator j=sp_nljm[norp].begin();j!=sp_nljm[norp].end();j++)
	// 		{
				
	// 			cout<<yq[distance(sp_nljm[norp].begin(),i)+distance(sp_nljm[norp].begin(),j)*sp_nljm[norp].size()]<<' ';
				
	// 		}
	// 		cout<<endl;
	// 	}
	// }
	PQ_M_struct temp;
	temp.J=P_Jpi[norp][P_i].J;
	temp.pi=P_Jpi[norp][P_i].pi;
	temp.PQ_Jpi_index=P_i;
	temp.M=M;
	temp.yq=yq;
	if(M!=0)
	{
		M=-M;
		yq=new double [sp_nljm_dim2[norp]];
		memset(yq,0,sizeof(double)*sp_nljm_dim2[norp]);
		for(vector<nljm_struct>::iterator i=sp_nljm[norp].begin();i!=sp_nljm[norp].end();i++)
		{
			for(vector<nljm_struct>::iterator j=sp_nljm[norp].begin();j!=sp_nljm[norp].end();j++)
			{
				if(i->m+j->m==M&&P_Jpi[norp][P_i].J>=abs(i->j-j->j)&&P_Jpi[norp][P_i].J<=i->j+j->j&&(i->l+j->l)%2==P_Jpi[norp][P_i].pi)
				{
					yq[distance(sp_nljm[norp].begin(),i)+distance(sp_nljm[norp].begin(),j)*sp_nljm[norp].size()]+=P_Jpi[norp][P_i].yq[i->nlj*sp_nlj[norp].size()+j->nlj]*cg_cal(i->j,i->m,j->j,j->m,P_Jpi[norp][P_i].J,M);
					//cout<<i->j<<' '<<i->m<<' '<<j->j<<' '<<j->m<<' '<<P_Jpi[norp][P_i].J<<' '<<M<<' '<<cg_cal(i->j,i->m,j->j,j->m,P_Jpi[norp][P_i].J,M)<<' '<<P_Jpi[norp][P_i].yq[i->nlj*sp_nlj[norp].size()+j->nlj]<<endl;
				}
			}
		}
		P_anti_symmetric(yq,norp);
		
		
		temp.yq_mM=yq;
	}
	else
	{
		temp.yq_mM=yq;
	}
	
	P_M[norp].emplace_back(temp);
	return P_M[norp].size()-1;
}

int QJpi_to_QM(int Q_i,int M,int norp)
{
	double *yq=new double [sp_nljm_dim2[norp]];
	memset(yq,0,sizeof(double)*sp_nljm_dim2[norp]);
	for(vector<nljm_struct>::iterator i=sp_nljm[norp].begin();i!=sp_nljm[norp].end();i++)
	{
		for(vector<nljm_struct>::iterator j=sp_nljm[norp].begin();j!=sp_nljm[norp].end();j++)
		{
			if(i->m-j->m==M&&Q_Jpi[norp][Q_i].J>=abs(i->j-j->j)&&Q_Jpi[norp][Q_i].J<=i->j+j->j&&(i->l+j->l)%2==Q_Jpi[norp][Q_i].pi)
			{
				if((j->j+j->m)%4==0)
				{
					yq[distance(sp_nljm[norp].begin(),i)+distance(sp_nljm[norp].begin(),j)*sp_nljm[norp].size()]+=Q_Jpi[norp][Q_i].yq[i->nlj*sp_nlj[norp].size()+j->nlj]*cg_cal(i->j,i->m,j->j,-j->m,Q_Jpi[norp][Q_i].J,M);
				}
				else
				{
					yq[distance(sp_nljm[norp].begin(),i)+distance(sp_nljm[norp].begin(),j)*sp_nljm[norp].size()]-=Q_Jpi[norp][Q_i].yq[i->nlj*sp_nlj[norp].size()+j->nlj]*cg_cal(i->j,i->m,j->j,-j->m,Q_Jpi[norp][Q_i].J,M);
				}
				
			}
		}
	}
	PQ_M_struct temp;
	temp.J=Q_Jpi[norp][Q_i].J;
	temp.pi=Q_Jpi[norp][Q_i].pi;
	temp.PQ_Jpi_index=Q_i;
	temp.M=M;
	temp.yq=yq;
	temp.yq_mM=new double [sp_nljm_dim2[norp]];
	memcpy(temp.yq_mM,yq,sizeof(double)*sp_nljm_dim2[norp]);
	mkl_dimatcopy('c','t',sp_nljm[norp].size(),sp_nljm[norp].size(),1,temp.yq_mM,sp_nljm[norp].size(),sp_nljm[norp].size());
	Q_M[norp].emplace_back(temp);
	return Q_M[norp].size()-1;
}
PQ_mat_struct PQ_mat_out(int J,int pi,int M,int PQ_Jpi_index,int PQ_M_index)
{
	PQ_mat_struct PQ;
	PQ.J=J;
	PQ.pi=pi;
	PQ.M=M;
	PQ.PQ_Jpi_index=PQ_Jpi_index;
	PQ.val=new int [2];
	PQ.val[0]=1;
	PQ.val[1]=1;
	PQ.ser=new int *[2];
	PQ.ser[0]=NULL;
	PQ.ser[1]=new int [3];
	PQ.ser[1][0]=2;
	PQ.ser[1][2]=1;
	PQ.ser[1][1]=PQ_M_index;
	return PQ;
}

PQ_mat_struct PQ_mat_out(int J,int pi,int M,int PQ_Jpi_index,int PQ_M_index,int lorr)
//this function is for the s.p. vector creation
{
	PQ_mat_struct PQ;
	PQ.J=J;
	PQ.pi=pi;
	PQ.M=M;
	PQ.PQ_Jpi_index=PQ_Jpi_index;
	PQ.val=new int [2];
	PQ.val[0]=1;
	PQ.val[1]=1;
	PQ.ser=new int *[2];
	PQ.ser[0]=NULL;
	PQ.ser[1]=new int [3];
	PQ.ser[1][0]=2;
	PQ.ser[1][2]=lorr;//lorr=0 the s.p. is a raw vector on the right side; ==1 s.p. a column vector on the left side;
	PQ.ser[1][1]=PQ_M_index;
	return PQ;
}
vector<vector<double *>> pow_gamma[2];
vector<vector<double *>> pow_mgamma[2];
vector<vector<double>> tr_pow_gamma[2];
vector<vector<double>> tr_pow_mgamma[2];
void pow_gamma_out(int norp)
{
	for(int i=0;i<P_M[norp].size();i++)
	{
		vector<double *> pow_gamma_temp;
		vector<double *> pow_mgamma_temp;
		vector<double> tr_pow_gamma_temp;
		vector<double> tr_pow_mgamma_temp;
		if(i==0)
		{
			pow_gamma[norp].emplace_back(pow_gamma_temp);
			pow_mgamma[norp].emplace_back(pow_mgamma_temp);
			tr_pow_gamma[norp].emplace_back(tr_pow_gamma_temp);
			tr_pow_mgamma[norp].emplace_back(tr_pow_mgamma_temp);
			continue;
		}
		
		for(int k=0;k<=N[norp]/2*2+2;k++)
		{
			double *pow_gamma_k_temp=new double [sp_nljm_dim2[norp]];
			double tr_pow_gamma_k_temp;
			double *y=P_M[norp][i].yq;
			if(k==0)
			{
				memset(pow_gamma_k_temp,0,sizeof(double)*sp_nljm_dim2[norp]);
				for(int p=0;p<sp_nljm[norp].size();p++)
				{
					pow_gamma_k_temp[p+p*sp_nljm[norp].size()]=1;
				}
				tr_pow_gamma_k_temp=sp_nljm[norp].size();		
			}
			else if(k==1)
			{
				memcpy(pow_gamma_k_temp,y,sizeof(double)*sp_nljm_dim2[norp]);
				tr_pow_gamma_k_temp=0;	
			}
			else if(k==2)
			{
				char ta='n';
				double alpha=1;
				double beta=0;
				int dim=sp_nljm[norp].size();
				dgemm_(&ta,&ta,&dim,&dim,&dim,&alpha,y,&dim,y,&dim,&beta,pow_gamma_k_temp,&dim);
				tr_pow_gamma_k_temp=0;
				for(int p=0;p<sp_nljm[norp].size();p++)
				{
					tr_pow_gamma_k_temp+=pow_gamma_k_temp[p+p*sp_nljm[norp].size()];
				}
			}
			else
			{
				char ta='n';
				double alpha=1;
				double beta=0;
				int dim=sp_nljm[norp].size();
				char side='l';
				char uorl='u';
				dsymm_(&side,&uorl,&dim,&dim,&alpha,pow_gamma_temp[2],&dim,pow_gamma_temp[k-2],&dim,&beta,pow_gamma_k_temp,&dim);
				tr_pow_gamma_k_temp=0;
				if(k%2==0)
				{
					for(int p=0;p<sp_nljm[norp].size();p++)
					{
						tr_pow_gamma_k_temp+=pow_gamma_k_temp[p+p*sp_nljm[norp].size()];
					}
				}
				
			}
			pow_gamma_temp.emplace_back(pow_gamma_k_temp);
			tr_pow_gamma_temp.emplace_back(tr_pow_gamma_k_temp);
		}
		pow_gamma[norp].emplace_back(pow_gamma_temp);
		tr_pow_gamma[norp].emplace_back(tr_pow_gamma_temp);
		if(P_M[norp][i].M!=0)
		{
			for(int k=0;k<=N[norp]/2*2+2;k++)
			{
				double *pow_mgamma_k_temp=new double [sp_nljm_dim2[norp]];
				double *y=P_M[norp][i].yq_mM;
				double tr_pow_mgamma_k_temp;
				if(k==0)
				{
					memset(pow_mgamma_k_temp,0,sizeof(double)*sp_nljm_dim2[norp]);
					for(int p=0;p<sp_nljm[norp].size();p++)
					{
						pow_mgamma_k_temp[p+p*sp_nljm[norp].size()]=1;
					}	
					tr_pow_mgamma_k_temp=sp_nljm[norp].size();
				}
				else if(k==1)
				{
					memcpy(pow_mgamma_k_temp,y,sizeof(double)*sp_nljm_dim2[norp]);
					tr_pow_mgamma_k_temp=0;
				}
				else if(k==2)
				{
					char ta='n';
					double alpha=1;
					double beta=0;
					int dim=sp_nljm[norp].size();
					dgemm_(&ta,&ta,&dim,&dim,&dim,&alpha,y,&dim,y,&dim,&beta,pow_mgamma_k_temp,&dim);
					tr_pow_mgamma_k_temp=0;
					for(int p=0;p<sp_nljm[norp].size();p++)
					{
						tr_pow_mgamma_k_temp+=pow_mgamma_k_temp[p+p*sp_nljm[norp].size()];
					}
				}
				else
				{
					double alpha=1;
					double beta=0;
					int dim=sp_nljm[norp].size();
					char side='l';
					char uorl='u';
					dsymm_(&side,&uorl,&dim,&dim,&alpha,pow_mgamma_temp[2],&dim,pow_mgamma_temp[k-2],&dim,&beta,pow_mgamma_k_temp,&dim);
					tr_pow_mgamma_k_temp=0;
					if(k%2==0)
					{
						for(int p=0;p<sp_nljm[norp].size();p++)
						{
							tr_pow_mgamma_k_temp+=pow_mgamma_k_temp[p+p*sp_nljm[norp].size()];
						}
					}
					
				}
				pow_mgamma_temp.emplace_back(pow_mgamma_k_temp);
				tr_pow_mgamma_temp.emplace_back(tr_pow_mgamma_k_temp);
			}
		}
		pow_mgamma[norp].emplace_back(pow_mgamma_temp);
		tr_pow_mgamma[norp].emplace_back(tr_pow_mgamma_temp);
	}
}

// size_t P_M_to_mat(int k,int norp)
// {
// 	PQ_mat_struct P_mat_temp;
// 	if(k>=0)
// 	{
// 		P_mat_temp.J=P_M[norp][k].J;
// 		P_mat_temp.M=P_M[norp][k].M;
// 		P_mat_temp.pi=P_M[norp][k].pi;
// 		vector<int> ser_temp;
// 		ser_temp.emplace_back(k);
// 		ser_temp.emplace_back(1);
// 		ser_int_struct mat_temp;
// 		mat_temp.val=1;
// 		mat_temp.ser=ser_temp;
// 		P_mat_temp.mat_tab.emplace_back(mat_temp);
// 		double coe;
// 		return PQ_mat_search_insert(P_mat[norp],P_mat_order[norp],P_mat_temp,coe);
// 	}
// 	else
// 	{
// 		P_mat_temp.J=P_M[norp][-k].J;
// 		P_mat_temp.M=-P_M[norp][-k].M;
// 		P_mat_temp.pi=P_M[norp][-k].pi;
// 		vector<int> ser_temp;
// 		ser_temp.emplace_back(k);
// 		ser_temp.emplace_back(1);
// 		ser_int_struct mat_temp;
// 		mat_temp.val=1;
// 		mat_temp.ser=ser_temp;
// 		P_mat_temp.mat_tab.emplace_back(mat_temp);
// 		double coe;
// 		return PQ_mat_search_insert(P_mat[norp],P_mat_order[norp],P_mat_temp,coe);
// 	}
// }






vector<vector<double *>> pow_q[2];
vector<vector<double *>> pow_qT[2];
vector<vector<double>> tr_pow_q[2];
void pow_q_out(int norp)
{
	for(int i=0;i<Q_M[norp].size();i++)
	{
		vector<double *> pow_q_temp;
        vector<double *> pow_qT_temp;
		vector<double> tr_pow_q_temp;
		if(i==0)
		{
			pow_q[norp].emplace_back(pow_q_temp);
            pow_qT[norp].emplace_back(pow_qT_temp);
			tr_pow_q[norp].emplace_back(tr_pow_q_temp);
			continue;
		}
		for(int k=0;k<=4;k++)
		{
			double *pow_q_k_temp=new double [sp_nljm_dim2[norp]];
			double tr_pow_q_k_temp;
            double *pow_qT_k_temp=new double [sp_nljm_dim2[norp]];
			double *y=Q_M[norp][i].yq;
            double *yT=Q_M[norp][i].yq_mM;
			if(k==0)
			{
				memset(pow_q_k_temp,0,sizeof(double)*sp_nljm_dim2[norp]);
                memset(pow_qT_k_temp,0,sizeof(double)*sp_nljm_dim2[norp]);
				for(int p=0;p<sp_nljm[norp].size();p++)
				{
					pow_q_k_temp[p+p*sp_nljm[norp].size()]=1;
                    pow_qT_k_temp[p+p*sp_nljm[norp].size()]=1;
				}
				tr_pow_q_k_temp=sp_nljm[norp].size();	
			}
			else if(k==1)
			{
				memcpy(pow_q_k_temp,y,sizeof(double)*sp_nljm_dim2[norp]);
                memcpy(pow_qT_k_temp,yT,sizeof(double)*sp_nljm_dim2[norp]);
				tr_pow_q_k_temp=0;
				for(int p=0;p<sp_nljm[norp].size();p++)
				{
					tr_pow_q_k_temp+=pow_q_k_temp[p+p*sp_nljm[norp].size()];
				}
			}
			else
			{
				char ta='n';
				double alpha=1;
				double beta=0;
				int dim=sp_nljm[norp].size();
				dgemm_(&ta,&ta,&dim,&dim,&dim,&alpha,pow_q_temp[k-1],&dim,y,&dim,&beta,pow_q_k_temp,&dim);
                dgemm_(&ta,&ta,&dim,&dim,&dim,&alpha,pow_qT_temp[k-1],&dim,yT,&dim,&beta,pow_qT_k_temp,&dim);
				tr_pow_q_k_temp=0;
				for(int p=0;p<sp_nljm[norp].size();p++)
				{
					tr_pow_q_k_temp+=pow_q_k_temp[p+p*sp_nljm[norp].size()];
				}
			}
			pow_q_temp.emplace_back(pow_q_k_temp);
            pow_qT_temp.emplace_back(pow_qT_k_temp);
			tr_pow_q_temp.emplace_back(tr_pow_q_k_temp);
		}
		pow_q[norp].emplace_back(pow_q_temp);	
        pow_qT[norp].emplace_back(pow_qT_temp);	
		tr_pow_q[norp].emplace_back(tr_pow_q_temp);
	}
}

void qT(int &val_in,int *key_in,int &val_out,int *key_out)
{
    val_out=val_in;
	int size=key_in[0];
	key_out[0]=size;
    if(size>2)
    {		
        for(int k=1;k<=size;k+=2)
        {
			key_out[k]=key_in[size-k];
			key_out[k+1]=key_in[size-k+1];
        }
	}
	else
	{
		memcpy(key_out+1,key_in+1,sizeof(int)*2);
	}
	
	for(int k=1;k<=size;k+=2)
	{
		if(key_out[k]>INT_MAX/2)
		{
			key_out[k]=-INT_MAX/2-(key_out[k]-INT_MAX/2);
		}
		else if(key_out[k]<-INT_MAX/2)
		{
			key_out[k]=INT_MAX/2-(key_out[k]+INT_MAX/2);
		}
		else if(abs(key_out[k])<INT_MAX/2&&abs(key_out[k])>INT_MAX/4)
		{
			key_out[k+1]=1-key_out[k+1];
		}
		else
		{
			if(key_out[k+1]%2!=0)
			{
				val_out*=-1;
			}
		}
	}        
}
void CT(PQ_mat_struct &C)
{
	int *val=new int [C.val[0]+1];
	int **ser=new int* [C.val[0]+1];
	val[0]=C.val[0];
	for(int i=1;i<=C.val[0];i++)
	{
		ser[i]=new int [C.ser[i][0]+1];
		qT(C.val[i],C.ser[i],val[i],ser[i]);
		delete [] C.ser[i];
		C.ser[i]=ser[i];
	}
	delete [] C.val;
	C.val=val;
	delete [] ser;
	//here C is freed and reevaluated, so don't worry about its mem leak.
}