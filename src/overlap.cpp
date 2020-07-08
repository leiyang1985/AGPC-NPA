#define N_check 2
double ove_cal(vector<PQ_mat_struct> &sr,int N_ove,int norp);
// bool operator < (const PQ_mat_struct &P1, const PQ_mat_struct &P2)
// {
// 	if(P1.M<P2.M)
// 	{
// 		return true;
// 	}
// 	else if(P1.M==P2.M)
// 	{
// 		if(P1.pi<P2.pi)
// 		{
// 			return true;
// 		}
// 		else if(P1.pi==P2.pi)
// 		{
// 			if(P1.mat_tab.size()<P2.mat_tab.size())
// 			{
// 				return true;
// 			}
// 			else if(P1.mat_tab.size()==P2.mat_tab.size())
// 			{
// 				for(int i=0;i<P1.mat_tab.size();i++)
// 				{
// 					if(P1.mat_tab[i].val<P1.mat_tab[i].val)
// 					{
// 						return true;
// 					}
// 					else if(P1.mat_tab[i].val==P2.mat_tab[i].val)
// 					{
// 						if(P1.mat_tab[i].ser.size()<P2.mat_tab[i].ser.size())
// 						{
// 							return true;
// 						}
// 						else if(P1.mat_tab[i].ser.size()==P2.mat_tab[i].ser.size())
// 						{
// 							if(memcmp(&(P1.mat_tab[i].ser[0]),&(P2.mat_tab[i].ser[0]),sizeof(int)*P1.mat_tab[i].ser.size())<0)
// 							{
// 								return true;
// 							}							
// 						}
// 					}
					
// 				}
// 			}
// 		}
// 	}
// 	return false;
// }


int PQ_mat_cmp(const void *p1,const void *p2)
{
	PQ_mat_struct *P1=(PQ_mat_struct *)p1;
	PQ_mat_struct *P2=(PQ_mat_struct *)p2;
	if(P1->M<P2->M)
	{
		return -1;
	}
	else if(P1->M>P2->M)
	{
		return 1;
	}
	else
	{
		if(P1->pi<P2->pi)
		{
			return -1;
		}
		else if(P1->pi>P2->pi)
		{
			return 1;
		}
		else
		{
			if(P1->val[0]<P2->val[0])
			{
				return -1;
			}
			else if(P1->val[0]>P2->val[0])
			{
				return 1;
			}
			else
			{
				for(int i=1;i<=P1->val[0];i++)
				{
					if(P1->val[i]<P2->val[i])
					{
						return -1;
					}
					else if(P1->val[i]>P2->val[i])
					{
						return 1;
					}
					else
					{
						if(P1->ser[i][0]<P2->ser[i][0])
						{
							return -1;
						}
						else if(P1->ser[i][0]>P2->ser[i][0])
						{
							return 1;
						}
						else
						{
							int test=memcmp(P1->ser[i]+1,P2->ser[i]+1,sizeof(int)*P1->ser[i][0]);
							if(test<0)
							{
								return -1;
							}	
							else if(test>0)
							{
								return 1;
							}						
						}
					}
					
				}
				return 0;
			}
		}
	}
}
void sr_sort(PQ_mat_struct * sr,int N_ove)
{
	qsort(sr,N_ove,sizeof(PQ_mat_struct),PQ_mat_cmp);
	qsort(sr+N_ove,N_ove,sizeof(PQ_mat_struct),PQ_mat_cmp);
	int diff_s=0;
	int diff_r=0;
	for(int i=0;i<N_ove-1;i++)
	{
		if(PQ_mat_cmp(sr+i,sr+i+1)!=0)
		{
			diff_s++;
		}
		if(PQ_mat_cmp(sr+i+N_ove,sr+i+1+N_ove)!=0)
		{
			diff_r++;
		}
	}
	if(diff_s>diff_r)
	{
		PQ_mat_struct temp[N_ove];
		memcpy(temp,sr,sizeof(PQ_mat_struct)*N_ove);
		memcpy(sr,sr+N_ove,sizeof(PQ_mat_struct)*N_ove);
		memcpy(sr+N_ove,temp,sizeof(PQ_mat_struct)*N_ove);
		return;
	}
	if(diff_s==diff_r)
	{
		for(int i=0;i<N_ove;i++)
		{
			if(PQ_mat_cmp(sr+i,sr+i+N_ove)>0)
			{
				PQ_mat_struct temp[N_ove];
				memcpy(temp,sr,sizeof(PQ_mat_struct)*N_ove);
				memcpy(sr,sr+N_ove,sizeof(PQ_mat_struct)*N_ove);
				memcpy(sr+N_ove,temp,sizeof(PQ_mat_struct)*N_ove);
				return;
			}
		}
	}
}
void sr_o_sort(PQ_mat_struct * sr,int N_ove)
{
	qsort(sr,N_ove,sizeof(PQ_mat_struct),PQ_mat_cmp);
	qsort(sr+N_ove+1,N_ove,sizeof(PQ_mat_struct),PQ_mat_cmp);
	int diff_s=0;
	int diff_r=0;
	for(int i=0;i<N_ove-1;i++)
	{
		if(PQ_mat_cmp(sr+i,sr+i+1)!=0)
		{
			diff_s++;
		}
		if(PQ_mat_cmp(sr+i+N_ove+1,sr+i+2+N_ove)!=0)
		{
			diff_r++;
		}
	}
	if(diff_s>diff_r)
	{
		PQ_mat_struct temp[N_ove+1];
		memcpy(temp,sr,sizeof(PQ_mat_struct)*(N_ove+1));
		memcpy(sr,sr+N_ove+1,sizeof(PQ_mat_struct)*(N_ove+1));
		memcpy(sr+N_ove+1,temp,sizeof(PQ_mat_struct)*(N_ove+1));
		CT(sr[N_ove]);
		CT(sr[2*N_ove+1]);
		return;
	}
	if(diff_s==diff_r)
	{
		for(int i=0;i<N_ove+1;i++)
		{
			if(PQ_mat_cmp(sr+i,sr+i+N_ove+1)>0)
			{
				PQ_mat_struct temp[N_ove+1];
				memcpy(temp,sr,sizeof(PQ_mat_struct)*(N_ove+1));
				memcpy(sr,sr+N_ove+1,sizeof(PQ_mat_struct)*(N_ove+1));
				memcpy(sr+N_ove+1,temp,sizeof(PQ_mat_struct)*(N_ove+1));
				CT(sr[N_ove]);
				CT(sr[2*N_ove+1]);
				return;
			}
		}
	}
}
int sr_cmp(PQ_mat_struct *&T1,PQ_mat_struct *&T2,int &N_ove)
{
	for(int i=0;i<2*N_ove;i++)
	{
		int test=PQ_mat_cmp(T1+i,T2+i);
		if(test<0)
		{
			return -1;
		}
		else if(test>0)
		{
			return 1;
		}
	}
	return 0;
}
int sr_o_cmp(PQ_mat_struct *&T1,PQ_mat_struct *&T2,int &N_ove)
{
	for(int i=0;i<2*N_ove+2;i++)
	{
		int test=PQ_mat_cmp(T1+i,T2+i);
		if(test<0)
		{
			return -1;
		}
		else if(test>0)
		{
			return 1;
		}
	}
	return 0;
}
template<typename T_key,typename T_val> double tree_ove_op(node<T_key,T_val> *&tree,int N_ove,int norp)
{
	// node<T_key,T_val> *p=tree;
	// stack<node<T_key,T_val>*>*stack=NULL;
	// double res=0;
	// while(p!=NULL||stack!=NULL)
	// {
	// 	while(p!=NULL)
	// 	{
	// 		push(stack,p);
	// 		p=p->l;
	// 	}
	// 	if(stack!=NULL)
	// 	{
	// 		p=top(stack);
	// 		if(fabs(p->v)>1e-13)
	// 		{
	// 			res+=p->v*ove_cal(p->k,N_ove,norp);
	// 		}
	// 		for(int i=0;i<2*N_ove;i++)
	// 		{
	// 			int size=p->k[i].val[0];
	// 			for(int k=1;k<=size;k++)
	// 			{
	// 				delete [] p->k[i].ser[k];
	// 			}
	// 			delete [] p->k[i].val;
	// 			delete [] p->k[i].ser;
	// 		}
	// 		delete [] p->k;
    //         free(p);
	// 		pop(stack);
	// 		p=p->r;
	// 	}
	// }
	// return res;
	
	if(tree!=NULL)
	{
		double res=tree_ove_op(tree->l,N_ove,norp);
		
		if(fabs(tree->v)>1e-13)
		{
			res+=tree->v*ove_cal(tree->k,N_ove,norp);
			// if(tree->k[0].ser[1][2]==1&&tree->k[1].ser[1][2]==5)
			// {
			// 	cin>>wat;
			// }
			// double v=tree->v;
			// double ove_v=ove_cal(tree->k,N_ove,norp);
			// res+=v*ove_v;
			// if(N_ove==6)
			// {
			// 	cout<<N_ove<<':';
			// 	for(int i=0;i<2*N_ove;i++)
			// 	{
			// 		int size=tree->k[i].val[0];
			// 		for(int k=1;k<=size;k++)
			// 		{
			// 			cout<<tree->k[i].ser[k][2]<<' ';
			// 		}
					
			// 	}
			// 	cout<<';'<<v<<' '<<ove_v<<endl;
			// }
		}
		for(int i=0;i<2*N_ove;i++)
		{
			int size=tree->k[i].val[0];
			for(int k=1;k<=size;k++)
			{
				delete [] tree->k[i].ser[k];
			}
			
			delete [] tree->k[i].val;
			delete [] tree->k[i].ser;
			// all++;
		}
		delete [] tree->k;
		
		res+=tree_ove_op(tree->r,N_ove,norp);
		free(tree);
		return res;
	}
	else
	{
		return 0;
	}
}

void ove_tree_search_add(struct node<PQ_mat_struct*,double> * &tree,PQ_mat_struct *k,double &v,int (*cmp)(PQ_mat_struct* &,PQ_mat_struct*&,int &),int N_ove)
{
	if(tree==NULL)
	{
		tree=(node<PQ_mat_struct*,double> *)malloc(sizeof(struct node<PQ_mat_struct*,double>));
		tree->l=NULL;
		tree->r=NULL;
		tree->k=new PQ_mat_struct [2*N_ove];
		memcpy(tree->k,k,sizeof(PQ_mat_struct)*2*N_ove);
        tree->v=v;
		
	}
    // else
    // {
    //     node<PQ_mat_struct *,double> *p=tree;
    //     node<PQ_mat_struct *,double> *p_last=NULL;
    //     bool isL=true;
    //     while(p!=NULL)
    //     {
    //         if(cmp(k,p->k,N_ove)==0)
    //         {
    //             p->v+=v;
    //             return;
    //         }
    //         else if(cmp(k,p->k,N_ove)<0)
    //         {
    //             p_last=p;
    //             p= p->l;
    //             isL=true;
    //         }
    //         else
    //         {
    //             p_last=p;
    //             p = p->r;
    //             isL=false;
    //         }
    //     }
    //     if(isL)
    //     {
    //         p_last->l=(node<PQ_mat_struct*,double> *)malloc(sizeof(struct node<PQ_mat_struct *,double>));
    //         p_last->l->l=NULL;
    //         p_last->l->r=NULL;
    //         p_last->l->k=new PQ_mat_struct [N_ove*2];
    //         memcpy(p_last->l->k,k,sizeof(PQ_mat_struct)*2*N_ove);
    //         p_last->l->v=v;
    //     }
    //     else
    //     {
    //         p_last->r=(node<PQ_mat_struct*,double> *)malloc(sizeof(struct node<PQ_mat_struct*,double>));
    //         p_last->r->l=NULL;
    //         p_last->r->r=NULL;
    //         p_last->r->k=new PQ_mat_struct [N_ove*2];
    //         memcpy(p_last->r->k,k,sizeof(PQ_mat_struct)*2*N_ove);
    //         p_last->r->v=v;
    //     }        
    // }
	else
    {
        if(cmp(k,tree->k,N_ove)==0)
        {
			for(int i=0;i<2*N_ove;i++)
			{
				int size=k[i].val[0];
				for(int p=1;p<=size;p++)
				{
					
					delete [] k[i].ser[p];
				}
				delete [] k[i].val;
				delete [] k[i].ser;
				// all++;
			}
            tree->v+=v;
        }
        else if(cmp(k,tree->k,N_ove)<0)
        {
            ove_tree_search_add(tree->l,k,v,cmp,N_ove);
        }
        else
        {
            ove_tree_search_add(tree->r,k,v,cmp,N_ove);
        }
    }
}
void ove_o_tree_search_add(struct node<PQ_mat_struct*,double> * &tree,PQ_mat_struct *k,double &v,int (*cmp)(PQ_mat_struct* &,PQ_mat_struct*&,int &),int N_ove)
{
	if(tree==NULL)
	{
		tree=(node<PQ_mat_struct*,double> *)malloc(sizeof(struct node<PQ_mat_struct*,double>));
		tree->l=NULL;
		tree->r=NULL;
		tree->k=new PQ_mat_struct [2*N_ove+2];
		memcpy(tree->k,k,sizeof(PQ_mat_struct)*(2*N_ove+2));
        tree->v=v;
		
	}
	else
    {
        if(cmp(k,tree->k,N_ove)==0)
        {
			for(int i=0;i<2*N_ove+2;i++)
			{
				int size=k[i].val[0];
				for(int p=1;p<=size;p++)
				{
					delete [] k[i].ser[p];
				}
				delete [] k[i].val;
				delete [] k[i].ser;
				// all++;
			}
            tree->v+=v;
        }
        else if(cmp(k,tree->k,N_ove)<0)
        {
            ove_o_tree_search_add(tree->l,k,v,cmp,N_ove);
        }
        else
        {
            ove_o_tree_search_add(tree->r,k,v,cmp,N_ove);
        }
    }
}
void PQ_mat_cpy(PQ_mat_struct *des,PQ_mat_struct *src,int size)
{
	for(int i=0;i<size;i++)
	{
		des[i].J=src[i].J;
		des[i].pi=src[i].pi;
		des[i].M=src[i].M;
		des[i].PQ_Jpi_index=src[i].PQ_Jpi_index;

		// des[i].val=new int [1];
		// des[i].ser=new int *[1];
		// des[i].val[0]=0;
		// des[i].ser[0]=NULL;
		//found++;
		
		des[i].val=new int [src[i].val[0]+1];
		memcpy(des[i].val,src[i].val,sizeof(int)*(src[i].val[0]+1));
		des[i].ser=new int *[src[i].val[0]+1];
		des[i].ser[0]=NULL;
		for(int k=1;k<=src[i].val[0];k++)
		{
			des[i].ser[k]=new int [src[i].ser[k][0]+1];
			memcpy(des[i].ser[k],src[i].ser[k],sizeof(int)*(src[i].ser[k][0]+1));
		}
	}
}
void PQ_mat_free(PQ_mat_struct &PQ)
{
	for(int k=1;k<=PQ.val[0];k++)
	{
		delete [] PQ.ser[k];
	}
	delete [] PQ.ser;
	delete [] PQ.val;
}
double ove_cal(PQ_mat_struct* sr,int N_ove,int norp)
//the user needs to keep the sr serial workable for overlap calculation (M and parity conservation) before the calling of the overlap function
{
	if (N_ove == 0)
	{
		return 1;
	}
	if (N_ove == 1)
	{
		return PP_dot_cal(sr[0], sr[1], norp);
		//  return 0;
	}
	PQ_mat_struct sr_low[(N_ove-1)*2];
	node<PQ_mat_struct*,double> *ove_tree=NULL;
	for (int i = 0; i < N_ove; i++)
	//int i=0;
	{
		// double phi =1;
		double phi= PP_dot_cal(sr[i], sr[N_ove * 2 - 1], norp);

		if(fabs(phi)>1e-13)
		{
			if (i == 0)
			{
				PQ_mat_cpy(sr_low, sr + 1,  (N_ove - 1) * 2);
				//copy(sr.begin()+1,sr.begin()+2*N_ove-1,sr_low.begin());
				// found+=(N_ove - 1) * 2;
			}
			else
			{
				PQ_mat_cpy(sr_low, sr,  i);
				PQ_mat_cpy(sr_low + i, sr + i + 1,  (N_ove - 1) * 2 - i);
				// found+=(N_ove - 1) * 2;
				//copy(sr.begin(),sr.begin()+i,sr_low.begin());
				//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
			}
			sr_sort(sr_low,N_ove-1);

			ove_tree_search_add(ove_tree,sr_low,phi,sr_cmp,N_ove-1);
		}		
	}
	int coe_pq;
	for (int i = 1; i < N_ove; i++)
	{
		auto Q = PP_con_cal(sr[i], sr[N_ove * 2 - 1]);
		// cout<<"new "<<all<<" Q"<<endl;
		// PQ_mat_struct Q;
		// Q.val=new int [1];
		// Q.val[0]=0;
		// Q.ser=new int *[1];
		// Q.ser[0]=NULL;
		for (int k = 0; k < i; k++)
		{
			if(k==0)
			{
				if(i==1)
				{
					PQ_mat_cpy(sr_low+1,sr+2,(N_ove-1)*2-1);
					//copy(sr.begin()+2,sr.begin()+2*N_ove-1,sr_low.begin()+1);
				}
				else
				{
					PQ_mat_cpy(sr_low + 1, sr+1,  i-1);
					PQ_mat_cpy(sr_low + i, sr + i + 1, (N_ove - 1) * 2 - i);
					//copy(sr.begin()+1,sr.begin()+i,sr_low.begin()+1);
					//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
				}
			}
			else
			{
				if(i==k+1)
				{
					PQ_mat_cpy(sr_low+1,sr,k);
					PQ_mat_cpy(sr_low+1+k,sr+k+2,(N_ove-1)*2-1-k);
					//copy(sr.begin(),sr.begin()+k,sr_low.begin()+1);
					//copy(sr.begin()+k+2,sr.begin()+2*N_ove-1,sr_low.begin()+k+1);
				}
				else
				{
					PQ_mat_cpy(sr_low+1,sr,k);
					PQ_mat_cpy(sr_low+1+k,sr+k+1,i-k-1);
					PQ_mat_cpy(sr_low+i,sr+i+1,(N_ove-1)*2-i);
					//copy(sr.begin(),sr.begin()+k,sr_low.begin()+1);
					//copy(sr.begin()+k+1,sr.begin()+i,sr_low.begin()+k+1);
					//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
				}
			}
			sr_low[0]=PQ_con_cal(sr[k], Q,coe_pq);
			// sr_low[0].val=new int [1];
			// sr_low[0].val[0]=0;
			// sr_low[0].ser=new int *[1];
			// sr_low[0].ser[0]=NULL;
			// if(N_ove==N_check)
			// {
			// 	cout<<sr_low[0]<<' '<<sr[0]<<' '<<PP_dot_out(sr_low[0],sr[0],0)<<endl;
			// 	cout<<sr_low[0]<<' '<<sr[1]<<' '<<PP_dot_out(sr_low[0],sr[1],0)<<endl;
			// }
			sr_sort(sr_low,N_ove-1);
			double coe=coe_pq;
			ove_tree_search_add(ove_tree,sr_low,coe,sr_cmp,N_ove-1);
		}

		for(int k=1;k<=Q.val[0];k++)
		{
			// found++;
			delete [] Q.ser[k];
		}
		// cout<<"delete "<<found<<" Q"<<endl;
		delete [] Q.ser;
		delete [] Q.val;
	}
	return tree_ove_op(ove_tree,N_ove-1,norp);
	//return 0;
}

vector<struct PQ_mat_struct *> *sr_ove[2];
vector<double> *phi_ove[2];
vector<double> *val_ove[2];
int *k_ove[2];
void ove_tail_init()
{
	for(int norp=0;norp<2;norp++)
	{
		sr_ove[norp]=new vector<struct PQ_mat_struct *> [N[norp]/2+1];
		phi_ove[norp]=new vector<double> [N[norp]/2+1];
		val_ove[norp]=new vector<double> [N[norp]/2+1];
		k_ove[norp]=new int [N[norp]/2+1];
		
		for(int i=0;i<=N[norp]/2;i++)
		{
			vector<struct PQ_mat_struct *> temp1;
			vector<double> temp2;
			PQ_mat_struct *temp3=new struct PQ_mat_struct [i*2];
			for(int k=0;k<i*2;k++)
			{
				temp3[k]=PQ_mat_out(0,0,0,-1,0);
			}
			
			temp1.emplace_back(temp3);
			temp2.emplace_back(0);
			sr_ove[norp][i]=temp1;
			phi_ove[norp][i]=temp2;
			val_ove[norp][i]=temp2;
			k_ove[norp][i]=0;
		}
	}
}
template<typename T_key,typename T_val> void tree_ove_to_sr_ove(node<T_key,T_val> *&tree,int N_ove,int norp)
{
	if(tree!=NULL)
	{
		tree_ove_to_sr_ove(tree->l,N_ove,norp);
		if(fabs(tree->v)>1e-13)
		{
			phi_ove[norp][N_ove].emplace_back(tree->v);
			sr_ove[norp][N_ove].emplace_back(tree->k);
			val_ove[norp][N_ove].emplace_back(0);
		}
		tree_ove_to_sr_ove(tree->r,N_ove,norp);
		free(tree);
	}
}
void ove_tail(int N_ove,int N_max,int norp);
double ove_cal_tail(PQ_mat_struct* sr,int N_ove,int norp)
{
	if(N_ove==0)
	{
		return 1;
	}
	else if(N_ove==1)
	{
		return PP_dot_cal(sr[0],sr[1],norp);
	}
	else
	{
		if(sr_ove[norp][N_ove].size()>0)
		{
			sr_ove[norp][N_ove].clear();
			val_ove[norp][N_ove].clear();
			phi_ove[norp][N_ove].clear();
		}
		sr_ove[norp][N_ove].emplace_back(sr);
		phi_ove[norp][N_ove].emplace_back(1);
		val_ove[norp][N_ove].emplace_back(0);
		k_ove[norp][N_ove]=0;
		ove_tail(N_ove,N_ove,norp);
		return val_ove[norp][N_ove][0];
	}
}
void ove_tail(int N_ove,int N_max,int norp)
//the user needs to keep the sr serial workable for overlap calculation (M and parity conservation) before the calling of the overlap function
{
	if (N_ove == 1)
	{
		val_ove[norp][2][k_ove[norp][2]]+=phi_ove[norp][1][k_ove[norp][N_ove]]*PP_dot_cal(sr_ove[norp][1][k_ove[norp][N_ove]][0],sr_ove[norp][1][k_ove[norp][N_ove]][1],norp);
		if(k_ove[norp][1]<phi_ove[norp][1].size()-1)
		{
			k_ove[norp][1]++;
			ove_tail(1,N_max,norp);
		}
		else
		{
			k_ove[norp][2]++;
			ove_tail(2,N_max,norp);
		}
	}
	else if(N_ove==N_max)
	{
		if(k_ove[norp][N_ove]>0)
		{
			return;
		}
		PQ_mat_struct sr_low[(N_ove-1)*2];
		node<PQ_mat_struct*,double> *ove_tree=NULL;
		for (int i = 0; i < N_ove; i++)
		//int i=0;
		{
			// double phi =1;
			double phi= PP_dot_cal(sr_ove[norp][N_ove][k_ove[norp][N_ove]][i], sr_ove[norp][N_ove][k_ove[norp][N_ove]][N_ove * 2 - 1], norp);
			if(fabs(phi)>1e-13)
			{
				if (i == 0)
				{
					PQ_mat_cpy(sr_low, sr_ove[norp][N_ove][k_ove[norp][N_ove]] + 1,  (N_ove - 1) * 2);
					//copy(sr.begin()+1,sr.begin()+2*N_ove-1,sr_low.begin());
					// found+=(N_ove - 1) * 2;
				}
				else
				{
					PQ_mat_cpy(sr_low, sr_ove[norp][N_ove][k_ove[norp][N_ove]],  i);
					PQ_mat_cpy(sr_low + i, sr_ove[norp][N_ove][k_ove[norp][N_ove]] + i + 1,  (N_ove - 1) * 2 - i);
					// found+=(N_ove - 1) * 2;
					//copy(sr.begin(),sr.begin()+i,sr_low.begin());
					//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
				}
				sr_sort(sr_low,N_ove-1);

				ove_tree_search_add(ove_tree,sr_low,phi,sr_cmp,N_ove-1);
			}		
		}
		int coe_pq;
		for (int i = 1; i < N_ove; i++)
		{
			auto Q = PP_con_cal(sr_ove[norp][N_ove][k_ove[norp][N_ove]][i], sr_ove[norp][N_ove][k_ove[norp][N_ove]][N_ove * 2 - 1]);
			// cout<<"new "<<all<<" Q"<<endl;
			// PQ_mat_struct Q;
			// Q.val=new int [1];
			// Q.val[0]=0;
			// Q.ser=new int *[1];
			// Q.ser[0]=NULL;
			for (int k = 0; k < i; k++)
			{
				if(k==0)
				{
					if(i==1)
					{
						PQ_mat_cpy(sr_low+1,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+2,(N_ove-1)*2-1);
						//copy(sr.begin()+2,sr.begin()+2*N_ove-1,sr_low.begin()+1);
					}
					else
					{
						PQ_mat_cpy(sr_low + 1, sr_ove[norp][N_ove][k_ove[norp][N_ove]]+1,  i-1);
						PQ_mat_cpy(sr_low + i, sr_ove[norp][N_ove][k_ove[norp][N_ove]] + i + 1, (N_ove - 1) * 2 - i);
						//copy(sr.begin()+1,sr.begin()+i,sr_low.begin()+1);
						//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
					}
				}
				else
				{
					if(i==k+1)
					{
						PQ_mat_cpy(sr_low+1,sr_ove[norp][N_ove][k_ove[norp][N_ove]],k);
						PQ_mat_cpy(sr_low+1+k,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+k+2,(N_ove-1)*2-1-k);
						//copy(sr.begin(),sr.begin()+k,sr_low.begin()+1);
						//copy(sr.begin()+k+2,sr.begin()+2*N_ove-1,sr_low.begin()+k+1);
					}
					else
					{
						PQ_mat_cpy(sr_low+1,sr_ove[norp][N_ove][k_ove[norp][N_ove]],k);
						PQ_mat_cpy(sr_low+1+k,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+k+1,i-k-1);
						PQ_mat_cpy(sr_low+i,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+i+1,(N_ove-1)*2-i);
						//copy(sr.begin(),sr.begin()+k,sr_low.begin()+1);
						//copy(sr.begin()+k+1,sr.begin()+i,sr_low.begin()+k+1);
						//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
					}
				}
				sr_low[0]=PQ_con_cal(sr_ove[norp][N_ove][k_ove[norp][N_ove]][k], Q,coe_pq);
				// sr_low[0].val=new int [1];
				// sr_low[0].val[0]=0;
				// sr_low[0].ser=new int *[1];
				// sr_low[0].ser[0]=NULL;
				// if(N_ove==N_check)
				// {
				// 	cout<<sr_low[0]<<' '<<sr[0]<<' '<<PP_dot_out(sr_low[0],sr[0],0)<<endl;
				// 	cout<<sr_low[0]<<' '<<sr[1]<<' '<<PP_dot_out(sr_low[0],sr[1],0)<<endl;
				// }
				sr_sort(sr_low,N_ove-1);
				double coe=coe_pq;
				ove_tree_search_add(ove_tree,sr_low,coe,sr_cmp,N_ove-1);
			}

			for(int k=1;k<=Q.val[0];k++)
			{
				// found++;
				delete [] Q.ser[k];
			}
			// cout<<"delete "<<found<<" Q"<<endl;
			delete [] Q.ser;
			delete [] Q.val;
		}
		for(int k=0;k<sr_ove[norp][N_ove-1].size();k++)
		{
			for(int i=0;i<2*(N_ove-1);i++)
			{
				PQ_mat_free(sr_ove[norp][N_ove-1][k][i]);
			}
			delete [] sr_ove[norp][N_ove-1][k];
		}
		sr_ove[norp][N_ove-1].clear();
		val_ove[norp][N_ove-1].clear();
		phi_ove[norp][N_ove-1].clear();
		tree_ove_to_sr_ove(ove_tree,N_ove-1,norp);
		k_ove[norp][N_ove-1]=0;
		ove_tail(N_ove-1,N_max,norp);
	}
	else
	{
		if(k_ove[norp][N_ove]>val_ove[norp][N_ove].size()-1)
		{
			for(int k=0;k<val_ove[norp][N_ove].size();k++)
			{
				val_ove[norp][N_ove+1][k_ove[norp][N_ove+1]]+=phi_ove[norp][N_ove][k]*val_ove[norp][N_ove][k];
				// if(N_ove==3)
				// {
				// 	cout<<phi_ove[norp][N_ove][k]<<' '<<val_ove[norp][N_ove][k]<<endl;
				// }
			}
			k_ove[norp][N_ove+1]++;
			ove_tail(N_ove+1,N_max,norp);
		}
		else
		{
			PQ_mat_struct sr_low[(N_ove-1)*2];
			node<PQ_mat_struct*,double> *ove_tree=NULL;
			for (int i = 0; i < N_ove; i++)
			//int i=0;
			{
				// double phi =1;
				double phi= PP_dot_cal(sr_ove[norp][N_ove][k_ove[norp][N_ove]][i], sr_ove[norp][N_ove][k_ove[norp][N_ove]][N_ove * 2 - 1], norp);
				if(fabs(phi)>1e-13)
				{
					if (i == 0)
					{
						PQ_mat_cpy(sr_low, sr_ove[norp][N_ove][k_ove[norp][N_ove]] + 1,  (N_ove - 1) * 2);
						//copy(sr.begin()+1,sr.begin()+2*N_ove-1,sr_low.begin());
						// found+=(N_ove - 1) * 2;
					}
					else
					{
						PQ_mat_cpy(sr_low, sr_ove[norp][N_ove][k_ove[norp][N_ove]],  i);
						PQ_mat_cpy(sr_low + i, sr_ove[norp][N_ove][k_ove[norp][N_ove]] + i + 1,  (N_ove - 1) * 2 - i);
						// found+=(N_ove - 1) * 2;
						//copy(sr.begin(),sr.begin()+i,sr_low.begin());
						//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
					}
					sr_sort(sr_low,N_ove-1);

					ove_tree_search_add(ove_tree,sr_low,phi,sr_cmp,N_ove-1);
				}		
			}
			int coe_pq;
			for (int i = 1; i < N_ove; i++)
			{
				auto Q = PP_con_cal(sr_ove[norp][N_ove][k_ove[norp][N_ove]][i], sr_ove[norp][N_ove][k_ove[norp][N_ove]][N_ove * 2 - 1]);
				// cout<<"new "<<all<<" Q"<<endl;
				// PQ_mat_struct Q;
				// Q.val=new int [1];
				// Q.val[0]=0;
				// Q.ser=new int *[1];
				// Q.ser[0]=NULL;
				for (int k = 0; k < i; k++)
				{
					if(k==0)
					{
						if(i==1)
						{
							PQ_mat_cpy(sr_low+1,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+2,(N_ove-1)*2-1);
							//copy(sr.begin()+2,sr.begin()+2*N_ove-1,sr_low.begin()+1);
						}
						else
						{
							PQ_mat_cpy(sr_low + 1, sr_ove[norp][N_ove][k_ove[norp][N_ove]]+1,  i-1);
							PQ_mat_cpy(sr_low + i, sr_ove[norp][N_ove][k_ove[norp][N_ove]] + i + 1, (N_ove - 1) * 2 - i);
							//copy(sr.begin()+1,sr.begin()+i,sr_low.begin()+1);
							//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
						}
					}
					else
					{
						if(i==k+1)
						{
							PQ_mat_cpy(sr_low+1,sr_ove[norp][N_ove][k_ove[norp][N_ove]],k);
							PQ_mat_cpy(sr_low+1+k,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+k+2,(N_ove-1)*2-1-k);
							//copy(sr.begin(),sr.begin()+k,sr_low.begin()+1);
							//copy(sr.begin()+k+2,sr.begin()+2*N_ove-1,sr_low.begin()+k+1);
						}
						else
						{
							PQ_mat_cpy(sr_low+1,sr_ove[norp][N_ove][k_ove[norp][N_ove]],k);
							PQ_mat_cpy(sr_low+1+k,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+k+1,i-k-1);
							PQ_mat_cpy(sr_low+i,sr_ove[norp][N_ove][k_ove[norp][N_ove]]+i+1,(N_ove-1)*2-i);
							//copy(sr.begin(),sr.begin()+k,sr_low.begin()+1);
							//copy(sr.begin()+k+1,sr.begin()+i,sr_low.begin()+k+1);
							//copy(sr.begin()+i+1,sr.begin()+2*N_ove-1,sr_low.begin()+i);
						}
					}
					sr_low[0]=PQ_con_cal(sr_ove[norp][N_ove][k_ove[norp][N_ove]][k], Q,coe_pq);
					// sr_low[0].val=new int [1];
					// sr_low[0].val[0]=0;
					// sr_low[0].ser=new int *[1];
					// sr_low[0].ser[0]=NULL;
					// if(N_ove==N_check)
					// {
					// 	cout<<sr_low[0]<<' '<<sr[0]<<' '<<PP_dot_out(sr_low[0],sr[0],0)<<endl;
					// 	cout<<sr_low[0]<<' '<<sr[1]<<' '<<PP_dot_out(sr_low[0],sr[1],0)<<endl;
					// }
					sr_sort(sr_low,N_ove-1);
					double coe=coe_pq;
					ove_tree_search_add(ove_tree,sr_low,coe,sr_cmp,N_ove-1);
				}

				for(int k=1;k<=Q.val[0];k++)
				{
					// found++;
					delete [] Q.ser[k];
				}
				// cout<<"delete "<<found<<" Q"<<endl;
				delete [] Q.ser;
				delete [] Q.val;
			}
			for(int k=0;k<sr_ove[norp][N_ove-1].size();k++)
			{
				for(int i=0;i<2*(N_ove-1);i++)
				{
					PQ_mat_free(sr_ove[norp][N_ove-1][k][i]);
				}
				delete [] sr_ove[norp][N_ove-1][k];
			}
			sr_ove[norp][N_ove-1].clear();
			val_ove[norp][N_ove-1].clear();
			phi_ove[norp][N_ove-1].clear();
			tree_ove_to_sr_ove(ove_tree,N_ove-1,norp);
			k_ove[norp][N_ove-1]=0;
			ove_tail(N_ove-1,N_max,norp);
		}
	}
}
double q_cal(PQ_mat_struct* sr,PQ_mat_struct &Q,int N_ove,int norp);
//in this code the animation s.p. operator is represented by raw vector; while cration by column vector
double ove_o_cal(PQ_mat_struct* sr,int norp)
{
	if(N[norp]/2==0)
	{
		return CC_dot_cal(sr[0],sr[1],norp);
	}
	else
	{
		PQ_mat_struct sr_low[(N[norp]/2)*2];
		double res=0;
		double phi=CC_dot_cal(sr[N[norp]/2],sr[(N[norp]/2)*2+1],norp);
		if(fabs(phi)>1e-13)
		{
			PQ_mat_cpy(sr_low,sr,N[norp]/2);
			PQ_mat_cpy(sr_low+N[norp]/2,sr+N[norp]/2+1,N[norp]/2);
			sr_sort(sr_low,N[norp]/2);
			res+=phi*ove_cal(sr_low,N[norp]/2,norp);
		}
		PQ_mat_cpy(sr_low,sr,N[norp]/2);
		PQ_mat_cpy(sr_low+N[norp]/2,sr+N[norp]/2+1,N[norp]/2);
		auto Q=CC_con_cal(sr[N[norp]/2],sr[(N[norp]/2)*2+1]);
		res+=q_cal(sr_low,Q,N[norp]/2,norp);
		PQ_mat_free(Q);
		return res;
	}
}
template<typename T_key,typename T_val> double tree_ove_o_op(node<T_key,T_val> *&tree,int N_ove,int norp)
{
	if(tree!=NULL)
	{
		double res=tree_ove_o_op(tree->l,N_ove,norp);
		
		if(fabs(tree->v)>1e-13)
		{
			res+=tree->v*ove_o_cal(tree->k,norp);
			// if(tree->k[0].ser[1][2]==1&&tree->k[1].ser[1][2]==5)
			// {
			// 	cin>>wat;
			// }
			// double v=tree->v;
			// double ove_v=ove_cal(tree->k,N_ove,norp);
			// res+=v*ove_v;
			// if(N_ove==6)
			// {
			// 	cout<<N_ove<<':';
			// 	for(int i=0;i<2*N_ove;i++)
			// 	{
			// 		int size=tree->k[i].val[0];
			// 		for(int k=1;k<=size;k++)
			// 		{
			// 			cout<<tree->k[i].ser[k][2]<<' ';
			// 		}
					
			// 	}
			// 	cout<<';'<<v<<' '<<ove_v<<endl;
			// }
		}
		for(int i=0;i<2*N_ove+2;i++)
		{
			int size=tree->k[i].val[0];
			for(int k=1;k<=size;k++)
			{
				delete [] tree->k[i].ser[k];
			}
			
			delete [] tree->k[i].val;
			delete [] tree->k[i].ser;
			// all++;
		}
		delete [] tree->k;
		
		res+=tree_ove_o_op(tree->r,N_ove,norp);
		free(tree);
		return res;
	}
	else
	{
		return 0;
	}
}