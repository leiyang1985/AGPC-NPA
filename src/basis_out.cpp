struct basis_block_struct
{
    int M;
    vector<vector<int>> r;
};
vector<basis_block_struct> basis_block[2][2];
//first is norp second is pi
int basis_M_max,basis_M_min;
int cmp_basis_block(basis_block_struct &T1,basis_block_struct &T2)
{
	if(T1.M<T2.M)
	{
		return -1;
	}
	else if(T1.M>T2.M)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}


void r_loop(int r[],int loop_index,int norp)
{
	int begin=0;
	if(loop_index>0)
	{
		begin=r[loop_index-1];
	}
	for(r[loop_index]=begin;r[loop_index]<P_in_basis[norp].size();r[loop_index]++)
	{
		if(loop_index==N[norp]/2-1)
		{
			// cout<<P_M[norp][abs(P_in_basis[norp][r[0]])].J<<' '<<P_M[norp][abs(P_in_basis[norp][r[0]])].M<<endl;
			for(int k=0;k<P_Jpi_in_basis_num[norp];k++)
			{
				if(P_num_max[norp][k]==-1)
				{
					continue;
				}
				int P_num=0;
				for(int i=0;i<N[norp]/2;i++)
				{
					if(P_M[norp][abs(P_in_basis[norp][r[i]])].PQ_Jpi_index==k)
					{
						P_num++;
						if(P_num>P_num_max[norp][k])
						{
							goto r_loop_next;
						}
					}
				}
			}
			if(N[norp]%2==0)
			{
				basis_block_struct basis_temp;
				basis_temp.M=0;
				int pi=0;
				for(int i=0;i<N[norp]/2;i++)
				{
					if(P_in_basis[norp][r[i]]>0)
					{
						basis_temp.M+=P_M[norp][P_in_basis[norp][r[i]]].M;
						pi+=P_M[norp][P_in_basis[norp][r[i]]].pi;
					}
					else
					{
						basis_temp.M-=P_M[norp][-P_in_basis[norp][r[i]]].M;
						pi+=P_M[norp][-P_in_basis[norp][r[i]]].pi;
					}
				}
				if(abs(basis_temp.M)>basis_M_max)
				{
					goto r_loop_next;
				}
				pi=pi%2;
				vector<int> r_insert;
				for(int i=0;i<N[norp]/2;i++)
				{
					r_insert.emplace_back(P_in_basis[norp][r[i]]);
				}
				qsort(&(r_insert[0]),N[norp]/2,sizeof(int),cmp);
				size_t index;

				if(array_search(&(basis_block[norp][pi][0]),basis_temp,index,basis_block[norp][pi].size(),cmp_basis_block))
				{
					size_t basis_index;
					array_search(basis_block[norp][pi][index].r,r_insert,basis_index,basis_block[norp][pi][index].r.size(),cmp);
					basis_block[norp][pi][index].r.emplace(basis_block[norp][pi][index].r.begin()+basis_index,r_insert);
				}
				else
				{
					// size_t basis_index;
					// array_search(basis_temp.r,r_insert,basis_index,basis_temp.r.size(),cmp);
					// basis_temp.r.emplace(basis_temp.r.begin()+basis_index,r_insert);
					basis_temp.r.emplace_back(r_insert);
					basis_block[norp][pi].emplace(basis_block[norp][pi].begin()+index,basis_temp);
				}
			}
			else
			{
				for(int k=0;k<sp_nljm[norp].size();k++)
				{
					if(sp_nlj[norp][sp_nljm[norp][k].nlj].unpair)
					{
						// cout<<sp_nljm[norp][k].j<<' '<<sp_nljm[norp][k].m<<endl;
						int sp_index;
						if(sp_nljm[norp][k].m>0)
						{
							sp_index=k+1;
						}
						else
						{
							sp_index=-(sp_nljm[norp].size()-k);
						}
						//the index of sp_nljm is 1 2 3 4 ...N -N -(N-1)... -4 -3 -2 -1 \pm INT_MAX/4 in the PQ_mat system
						basis_block_struct basis_temp;
						basis_temp.M=0;
						int pi=0;
						for(int i=0;i<N[norp]/2;i++)
						{
							if(P_in_basis[norp][r[i]]>0)
							{
								basis_temp.M+=P_M[norp][P_in_basis[norp][r[i]]].M;
								pi+=P_M[norp][P_in_basis[norp][r[i]]].pi;
							}
							else
							{
								basis_temp.M-=P_M[norp][-P_in_basis[norp][r[i]]].M;
								pi+=P_M[norp][-P_in_basis[norp][r[i]]].pi;
							}
						}
						if(sp_index>0)
						{
							basis_temp.M+=sp_nljm[norp][sp_index-1].m;
						}
						else
						{
							basis_temp.M-=sp_nljm[norp][-sp_index-1].m;
						}
						pi+=sp_nljm[norp][abs(sp_index)-1].l;
						if(abs(basis_temp.M)>basis_M_max)
						{
							continue;
						}
						pi=pi%2;
						vector<int> r_insert;
						for(int i=0;i<N[norp]/2;i++)
						{
							r_insert.emplace_back(P_in_basis[norp][r[i]]);
						}
						r_insert.emplace_back(sp_index);
						qsort(&(r_insert[0]),N[norp]/2,sizeof(int),cmp);
						size_t block_index;

						if(array_search(&(basis_block[norp][pi][0]),basis_temp,block_index,basis_block[norp][pi].size(),cmp_basis_block))
						{
							size_t basis_index;
							array_search(basis_block[norp][pi][block_index].r,r_insert,basis_index,basis_block[norp][pi][block_index].r.size(),cmp);
							basis_block[norp][pi][block_index].r.emplace(basis_block[norp][pi][block_index].r.begin()+basis_index,r_insert);
						}
						else
						{
							basis_temp.r.emplace_back(r_insert);
							basis_block[norp][pi].emplace(basis_block[norp][pi].begin()+block_index,basis_temp);
						}
					}
				}
			}
			
		}
		else
		{
			r_loop(r,loop_index+1,norp);
		}
r_loop_next: continue;
	}
}
int * Rev_p_M0[2][2];
void basis_out(int norp)
{
	int M_max[2]={0,0};
	for(int i=0;i<P_in_basis[0].size();i++)
	{
		int M=P_M[0][abs(P_in_basis[0][i])].M*(N[0]/2);
		if(N[0]%2==1)
		{
			int M_sp_max=0;
			for(int k=0;k<sp_nlj[0].size();k++)
			{
				if(sp_nlj[0][k].j>M_sp_max&&sp_nlj[0][k].unpair)
				{
					M_sp_max=sp_nlj[0][k].j;
				}
			}
			M+=M_sp_max;
		}
		if(M>M_max[0])
		{
			M_max[0]=M;
		}
	}
	for(int i=0;i<P_in_basis[1].size();i++)
	{
		int M=P_M[1][abs(P_in_basis[1][i])].M*(N[1]/2);
		if(N[1]%2==1)
		{
			int M_sp_max=0;
			for(int k=0;k<sp_nlj[1].size();k++)
			{
				if(sp_nlj[1][k].j>M_sp_max&&sp_nlj[1][k].unpair)
				{
					M_sp_max=sp_nlj[1][k].j;
				}
			}
			M+=M_sp_max;
		}
		if(M>M_max[1])
		{
			M_max[1]=M;
		}
	}
	basis_M_max=min(M_max[norp],M_max[1-norp]+M_tot);
	basis_M_min=max(-M_max[norp],-M_max[1-norp]+M_tot);
	if(abs(basis_M_min)>abs(basis_M_max))
	{
		basis_M_max=abs(basis_M_min);
	}
	else
	{
		basis_M_max=abs(basis_M_max);
	}
	//basis_M_max=20;
	if(N[norp]>1)
	{
		int r[N[norp]/2+N[norp]%2];
		r_loop(r,0,norp);
	}
	else if(N[norp]==1)
	{
		for(int k=0;k<sp_nljm[norp].size();k++)
		{
			if(sp_nlj[norp][sp_nljm[norp][k].nlj].unpair)
			{
				int sp_index;
				if(sp_nljm[norp][k].m>0)
				{
					sp_index=k+1;
				}
				else
				{
					sp_index=-(sp_nljm[norp].size()-k);
				}
				//the index of sp_nljm is 1 2 3 4 ...N -N -(N-1)... -4 -3 -2 -1 \pm INT_MAX/4 in the PQ_mat system
				basis_block_struct basis_temp;
				basis_temp.M=0;
				if(sp_index>0)
				{
					basis_temp.M=sp_nljm[norp][sp_index-1].m;
				}
				else
				{
					basis_temp.M=-sp_nljm[norp][-sp_index-1].m;
				}
				if(abs(basis_temp.M)>basis_M_max)
				{
					continue;
				}
				int pi=sp_nljm[norp][abs(sp_index)-1].l%2;
				vector<int> r_insert;
				r_insert.emplace_back(sp_index);
				size_t block_index;
				if(array_search(&(basis_block[norp][pi][0]),basis_temp,block_index,basis_block[norp][pi].size(),cmp_basis_block))
				{
					size_t basis_index;
					array_search(basis_block[norp][pi][block_index].r,r_insert,basis_index,basis_block[norp][pi][block_index].r.size(),cmp);
					basis_block[norp][pi][block_index].r.emplace(basis_block[norp][pi][block_index].r.begin()+basis_index,r_insert);
				}
				else
				{
					basis_temp.r.emplace_back(r_insert);
					basis_block[norp][pi].emplace(basis_block[norp][pi].begin()+block_index,basis_temp);
				}
			}
		}
	}
	else
	{
		basis_block_struct basis_temp;
		basis_temp.M=0;
		vector<int> r_insert;
		basis_temp.r.emplace_back(r_insert);
		basis_block[norp][0].emplace_back(basis_temp);
	}
	for(int pi=0;pi<2;pi++)
	{
		vector<int> s_mM(N[norp]/2+N[norp]%2);
		for(int i=0;i<basis_block[norp][pi].size()/2;i++)
		{
			int Rev_i=basis_block[norp][pi].size()-1-i;
			for(int p=0;p<basis_block[norp][pi][i].r.size();p++)
			{
				copy(basis_block[norp][pi][i].r[p].begin(),basis_block[norp][pi][i].r[p].end(),s_mM.begin());
				for(int k=0;k<N[norp]/2;k++)
				{
					if(P_M[norp][abs(s_mM[k])].M!=0)
					{
						s_mM[k]=-s_mM[k];
					}
				}
				if(N[norp]%2!=0)
				{
					s_mM[N[norp]/2]=-s_mM[N[norp]/2];
				}
				qsort(&(s_mM[0]),N[norp]/2,sizeof(int),cmp);
				copy(s_mM.begin(),s_mM.end(),basis_block[norp][pi][Rev_i].r[p].begin());
			}
		}
		if(basis_block[norp][pi].size()%2==1)
		{
			int i_M0=basis_block[norp][pi].size()/2;
			int dim_M0=basis_block[norp][pi][i_M0].r.size();
			Rev_p_M0[norp][pi]=new int [dim_M0];
			for(int p=0;p<dim_M0;p++)
			{
				copy(basis_block[norp][pi][i_M0].r[p].begin(),basis_block[norp][pi][i_M0].r[p].end(),s_mM.begin());
				for(int k=0;k<N[norp]/2;k++)
				{
					if(P_M[norp][abs(s_mM[k])].M!=0)
					{
						s_mM[k]=-s_mM[k];
					}
				}
				qsort(&(s_mM[0]),N[norp]/2,sizeof(int),cmp);
				size_t mM_index;
				array_search(basis_block[norp][pi][i_M0].r,s_mM,mM_index,basis_block[norp][pi][i_M0].r.size(),cmp);
				Rev_p_M0[norp][pi][p]=mM_index;
			}
		}
	}
}
