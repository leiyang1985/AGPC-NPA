struct ham_q_block_struct
{
    double ***Q_np[2];
    double ***Q_np_unnor[2];
    bool **Q_np_evalued[2];
    double *Jpm[2];
    double *Jpm_unnor[2];
    bool Jpm_evalued[2];
};
ham_q_block_struct **ham_q_block[2][2][2];
void mat_ort_nor(double *mat_out,double *mat, double *vec_l,double *vec_r, int dim_l,int dim_r)
{
	int dim2 = dim_l*dim_r;
	memset(mat_out, 0, sizeof(double) * dim2);
	for (int p = 0; p < dim_l; p++)
	{
		for (int q = 0; q < dim_r; q++)
		{
			for (int ii = 0; ii < p+1; ii++)
			{
				for (int jj = 0; jj < q+1; jj++)
				{
					mat_out[p+q*dim_l] += mat[ii+jj*dim_l] * vec_l[pack_k(ii, p)] * vec_r[pack_k(jj, q)];
				}
			}
		}
	}
}
void ham_q_block_insert(int i,int j,int k,int qi,int morp,int norp,int pi_l,int pi_r)
//morp=0: original operator; morp=1: hermiteed operator
{
    // if(i==12&&j==12&&k==1&&qi==0&&morp==0&&norp==0&&pi_l==0&&pi_r==0)
    // {
    //     cin>>wat;
    // }
    int her=QQ_np_ham[k].her;
    int Q_index=QQ_np_ham[k].Q_index[norp][qi];
    int M_l=basis_block[norp][pi_l][i].M;
    int M_r=basis_block[norp][pi_r][j].M;
    int dim_l=basis_block[norp][pi_l][i].r.size();
    int dim_r=basis_block[norp][pi_r][j].r.size();
    int Q_k=Q_M[norp][Q_index-INT_MAX/2].J;
    double *val;
    if(ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp][k][qi])
    {
        return;
    }
    
    int Rev_i=basis_block[norp][pi_l].size()-1-i;
    int Rev_j=basis_block[norp][pi_r].size()-1-j;
    if(ham_q_block[norp][pi_l][pi_r][Rev_i][Rev_j].Q_np_evalued[1-morp][k][qi])
    {
        val=new double [dim_l*dim_r];
        for(int p=0;p<basis_block[norp][pi_l][i].r.size();p++)
        {
            int Rev_p;
            if(basis_block[norp][pi_l][i].M!=0)
            {
                Rev_p=p;
            }
            else
            {
                Rev_p=Rev_p_M0[norp][pi_l][p];
            }
            int phase_l=0;
            for(int kk=0;kk<N[norp]/2;kk++)
            {
                phase_l+=P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J;
            }
            if(N[norp]%2==1)
            {
                phase_l+=sp_nljm[norp][abs(basis_block[norp][pi_l][i].r[p][N[norp]/2])-1].j;
            }
            phase_l-=M_l;
            for(int q=0;q<basis_block[norp][pi_r][j].r.size();q++)
            {
                int Rev_q;
                if(basis_block[norp][pi_r][j].M!=0)
                {
                    Rev_q=q;
                }
                else
                {
                    Rev_q=Rev_p_M0[norp][pi_r][q];
                }
                int phase=0;
                int phase_r=0;
                for(int kk=0;kk<N[norp]/2;kk++)
                {
                    phase_r+=P_M[norp][abs(basis_block[norp][pi_r][j].r[q][kk])].J;
                }
                if(N[norp]%2==1)
                {
                    phase_r+=sp_nljm[norp][abs(basis_block[norp][pi_r][j].r[q][N[norp]/2])-1].j;
                }
                phase_r-=M_r;
                phase=phase_l+phase_r+Q_k;
                if(phase%4==0)
                {
                    phase=1;
                }
                else
                {
                    phase=-1;
                }
                phase*=her;
                val[p+q*dim_l]=phase*ham_q_block[norp][pi_l][pi_r][Rev_i][Rev_j].Q_np_unnor[1-morp][k][qi][Rev_p+Rev_q*dim_l];
                // if(isnan(val[p+q*dim_l])==1)
                // {
                //     cin>>wat;
                // }
            }
        }
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi]=new double [dim_r*dim_l];
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[morp][k][qi]=val;
        mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp][k][qi]=true;
        return;
    }
    if(ham_q_block[norp][pi_r][pi_l][j][i].Q_np_evalued[1-morp][k][qi])
    {
        val=new double [dim_r*dim_l];
        memcpy(val,ham_q_block[norp][pi_r][pi_l][j][i].Q_np_unnor[1-morp][k][qi],sizeof(double)*dim_l*dim_r);
        mkl_dimatcopy('c', 't', dim_r, dim_l, 1.0, val, dim_r, dim_l);
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi]=new double [dim_r*dim_l];
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[morp][k][qi]=val;
        mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp][k][qi]=true;
        return;
    }
    
    if(ham_q_block[norp][pi_r][pi_l][Rev_j][Rev_i].Q_np_evalued[morp][k][qi])
    {
        val=new double [dim_l*dim_r];
        for(int p=0;p<basis_block[norp][pi_l][i].r.size();p++)
        {
            int Rev_p;
            if(basis_block[norp][pi_l][i].M!=0)
            {
                Rev_p=p;
            }
            else
            {
                Rev_p=Rev_p_M0[norp][pi_l][p];
            }
            int phase_l=0;
            for(int kk=0;kk<N[norp]/2;kk++)
            {
                phase_l+=P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J;
            }
            if(N[norp]%2==1)
            {
                phase_l+=sp_nljm[norp][abs(basis_block[norp][pi_l][i].r[p][N[norp]/2])-1].j;
            }
            phase_l-=M_l;
            for(int q=0;q<basis_block[norp][pi_r][j].r.size();q++)
            {
                int Rev_q;
                if(basis_block[norp][pi_r][j].M!=0)
                {
                    Rev_q=q;
                }
                else
                {
                    Rev_q=Rev_p_M0[norp][pi_r][q];
                }
                int phase=0;
                int phase_r=0;
                for(int kk=0;kk<N[norp]/2;kk++)
                {
                    phase_r+=P_M[norp][abs(basis_block[norp][pi_r][j].r[q][kk])].J;
                }
                if(N[norp]%2==1)
                {
                    phase_r+=sp_nljm[norp][abs(basis_block[norp][pi_r][j].r[q][N[norp]/2])-1].j;
                }
                phase_r-=M_r;
                phase=phase_l+phase_r+Q_k;
                if(phase%4==0)
                {
                    phase=1;
                }
                else
                {
                    phase=-1;
                }
                phase*=her;
                val[p+q*dim_l]=phase*ham_q_block[norp][pi_r][pi_l][Rev_j][Rev_i].Q_np_unnor[morp][k][qi][Rev_q+Rev_p*dim_r];
            }
        }
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi]=new double [dim_r*dim_l];
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[morp][k][qi]=val;
        mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp][k][qi]=true;
        return;
    }
    
    cout<<"block "<<i<<" in "<<basis_block[norp][pi_l].size()<<" of pi_l= "<<pi_l<<" ; "<<j<<" in "<<basis_block[norp][pi_r].size()<<" of pi_r= "<<pi_r<<" "<<norp<<"th neuclon k="<<k<<" qi="<<qi<<" morp="<<morp<<endl;
    PQ_mat_struct sr[N[norp]/2*2+N[norp]%2];
    val=new double [dim_r*dim_l];
    for(int p=0;p<dim_l;p++)
    {
         for(int kk=0;kk<N[norp]/2;kk++)
        {
            if(basis_block[norp][pi_l][i].r[p][kk]<0)
            {
                sr[kk]=PQ_mat_out(P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].pi,-P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].M,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].PQ_Jpi_index,basis_block[norp][pi_l][i].r[p][kk]);
            }
            else
            {
                sr[kk]=PQ_mat_out(P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].pi,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].M,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].PQ_Jpi_index,basis_block[norp][pi_l][i].r[p][kk]);
            }
        }
        if(N[norp]%2==1)
        {
            if(basis_block[norp][pi_l][i].r[p][N[norp]/2]>0)
            {
                int sp_index=basis_block[norp][pi_l][i].r[p][N[norp]/2]-1;
                sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,0);
            }
            else
            {
                int sp_index=-basis_block[norp][pi_l][i].r[p][N[norp]/2]-1;
                sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),0);
            }
        }
        for(int q=0;q<dim_r;q++)
        {
            for(int kk=0;kk<N[norp]/2;kk++)
            {
                if(basis_block[norp][pi_r][j].r[q][kk]<0)
                {
                    sr[kk+N[norp]/2+N[norp]%2]=PQ_mat_out(P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].J,P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].pi,-P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].M,P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].PQ_Jpi_index,basis_block[norp][pi_r][j].r[q][kk]);
                }
                else
                {
                    sr[kk+N[norp]/2+N[norp]%2]=PQ_mat_out(P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].J,P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].pi,P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].M,P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].PQ_Jpi_index,basis_block[norp][pi_r][j].r[q][kk]);
                }
            }
            if(N[norp]%2==1)
            {
                if(basis_block[norp][pi_r][j].r[q][N[norp]/2]>0)
                {
                    int sp_index=basis_block[norp][pi_r][j].r[q][N[norp]/2]-1;
                    sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,1);
                }
                else
                {
                    int sp_index=-basis_block[norp][pi_r][j].r[q][N[norp]/2]-1;
                    sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),1);
                }
            }
            PQ_mat_struct Q;
            if(morp==0)
            {
                // for(int i=0;i<sp_nljm[norp].size();i++)
                // {
                //     for(int j=0;j<sp_nljm[norp].size();j++)
                //     {
                //         cout<<Q_M[norp][Q_index-INT_MAX/2].yq[i+j*sp_nljm[norp].size()]<<' ';
                //     }
                //     cout<<endl;
                // }
                Q=PQ_mat_out(Q_M[norp][Q_index-INT_MAX/2].J,Q_M[norp][Q_index-INT_MAX/2].pi,Q_M[norp][Q_index-INT_MAX/2].M,Q_M[norp][Q_index-INT_MAX/2].PQ_Jpi_index,Q_index);
            }
            else
            {
                // for(int i=0;i<sp_nljm[norp].size();i++)
                // {
                //     for(int j=0;j<sp_nljm[norp].size();j++)
                //     {
                //         cout<<Q_M[norp][Q_index-INT_MAX/2].yq[i+j*sp_nljm[norp].size()]<<' ';
                //     }
                //     cout<<endl;
                // }
                Q=PQ_mat_out(Q_M[norp][Q_index-INT_MAX/2].J,Q_M[norp][Q_index-INT_MAX/2].pi,-Q_M[norp][Q_index-INT_MAX/2].M,Q_M[norp][Q_index-INT_MAX/2].PQ_Jpi_index,-Q_index);
            }
            if(N[norp]%2==0)
            {
                val[p+q*dim_l]=q_cal(sr,Q,N[norp]/2,norp);
                // if(isnan(val[p+q*dim_l])==1)
                // {
                //     cin>>wat;
                // }
            }
            else
            {
                val[p+q*dim_l]=q_o_cal(sr,Q,N[norp]/2,norp);
            }
            
            // cout<<val[p+q*dim_l]<<endl;
            PQ_mat_free(Q);
        }
        // cout<<endl;
    }
    ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi]=new double [dim_r*dim_l];
    ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[morp][k][qi]=val;
    mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k][qi],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
    ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp][k][qi]=true;
}

void Jpm_q_block_insert(int i,int j,int morp,int norp,int pi_l,int pi_r)
{
    if(ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[morp])
    {
        return;
    }
    int M_l=basis_block[norp][pi_l][i].M;
    int M_r=basis_block[norp][pi_r][j].M;
    int dim_l=basis_block[norp][pi_l][i].r.size();
    int dim_r=basis_block[norp][pi_r][j].r.size();
    double *val;
    
    int Rev_i=basis_block[norp][pi_l].size()-1-i;
    int Rev_j=basis_block[norp][pi_r].size()-1-j;
    if(ham_q_block[norp][pi_l][pi_r][Rev_i][Rev_j].Jpm_evalued[1-morp])
    {
        val=new double [dim_l*dim_r];
        for(int p=0;p<basis_block[norp][pi_l][i].r.size();p++)
        {
            int Rev_p;
            if(basis_block[norp][pi_l][i].M!=0)
            {
                Rev_p=p;
            }
            else
            {
                Rev_p=Rev_p_M0[norp][pi_l][p];
            }
            int phase_l=0;
            for(int kk=0;kk<N[norp]/2;kk++)
            {
                phase_l+=P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J;
            }
            if(N[norp]%2==1)
            {
                phase_l+=sp_nljm[norp][abs(basis_block[norp][pi_l][i].r[p][N[norp]/2])-1].j;
            }
            phase_l-=M_l;
            for(int q=0;q<basis_block[norp][pi_r][j].r.size();q++)
            {
                int Rev_q;
                if(basis_block[norp][pi_r][j].M!=0)
                {
                    Rev_q=q;
                }
                else
                {
                    Rev_q=Rev_p_M0[norp][pi_r][q];
                }
                int phase=0;
                int phase_r=0;
                for(int kk=0;kk<N[norp]/2;kk++)
                {
                    phase_r+=P_M[norp][abs(basis_block[norp][pi_r][j].r[q][kk])].J;
                }
                if(N[norp]%2==1)
                {
                    phase_r+=sp_nljm[norp][abs(basis_block[norp][pi_r][j].r[q][N[norp]/2])-1].j;
                }
                phase_r-=M_r;
                phase=phase_l+phase_r+2;
                if(phase%4==0)
                {
                    phase=1;
                }
                else
                {
                    phase=-1;
                }
                val[p+q*dim_l]=phase*ham_q_block[norp][pi_l][pi_r][Rev_i][Rev_j].Jpm_unnor[1-morp][Rev_p+Rev_q*dim_l];
            }
        }
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp]=new double [dim_r*dim_l];
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm_unnor[morp]=val;
        mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[morp]=true;
        return;
    }
    if(ham_q_block[norp][pi_r][pi_l][j][i].Jpm_evalued[1-morp])
    {
        val=new double [dim_l*dim_r];
        memcpy(val,ham_q_block[norp][pi_r][pi_l][j][i].Jpm_unnor[1-morp],sizeof(double)*dim_l*dim_r);
        mkl_dimatcopy('c', 't', dim_r, dim_l, 1.0, val, dim_r, dim_l);
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp]=new double [dim_r*dim_l];
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm_unnor[morp]=val;
        mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[morp]=true;
        return;
    }
    if(ham_q_block[norp][pi_r][pi_l][Rev_j][Rev_i].Jpm_evalued[morp])
    {
        val=new double [dim_l*dim_r];
        for(int p=0;p<basis_block[norp][pi_l][i].r.size();p++)
        {
            int Rev_p;
            if(basis_block[norp][pi_l][i].M!=0)
            {
                Rev_p=p;
            }
            else
            {
                Rev_p=Rev_p_M0[norp][pi_l][p];
            }
            int phase_l=0;
            for(int kk=0;kk<N[norp]/2;kk++)
            {
                phase_l+=P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J;
            }
            if(N[norp]%2==1)
            {
                phase_l+=sp_nljm[norp][abs(basis_block[norp][pi_l][i].r[p][N[norp]/2])-1].j;
            }
            phase_l-=M_l;
            for(int q=0;q<basis_block[norp][pi_r][j].r.size();q++)
            {
                int Rev_q;
                if(basis_block[norp][pi_r][j].M!=0)
                {
                    Rev_q=q;
                }
                else
                {
                    Rev_q=Rev_p_M0[norp][pi_r][q];
                }
                int phase=0;
                int phase_r=0;
                for(int kk=0;kk<N[norp]/2;kk++)
                {
                    phase_r+=P_M[norp][abs(basis_block[norp][pi_r][j].r[q][kk])].J;
                }
                if(N[norp]%2==1)
                {
                    phase_r+=sp_nljm[norp][abs(basis_block[norp][pi_r][j].r[q][N[norp]/2])-1].j;
                }
                phase_r-=M_r;
                phase=phase_l+phase_r+2;
                if(phase%4==0)
                {
                    phase=1;
                }
                else
                {
                    phase=-1;
                }
                val[p+q*dim_l]=phase*ham_q_block[norp][pi_r][pi_l][Rev_j][Rev_i].Jpm_unnor[morp][Rev_q+Rev_p*dim_r];
            }
        }
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp]=new double [dim_r*dim_l];
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm_unnor[morp]=val;
        mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[morp]=true;
        return;
    }
    
    PQ_mat_struct sr[N[norp]/2*2+N[norp]%2];
    val=new double [dim_r*dim_l];
    for(int p=0;p<dim_l;p++)
    {
         for(int kk=0;kk<N[norp]/2;kk++)
        {
            if(basis_block[norp][pi_l][i].r[p][kk]<0)
            {
                sr[kk]=PQ_mat_out(P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].pi,-P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].M,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].PQ_Jpi_index,basis_block[norp][pi_l][i].r[p][kk]);
            }
            else
            {
                sr[kk]=PQ_mat_out(P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].J,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].pi,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].M,P_M[norp][abs(basis_block[norp][pi_l][i].r[p][kk])].PQ_Jpi_index,basis_block[norp][pi_l][i].r[p][kk]);
            }
        }
        if(N[norp]%2==1)
        {
            if(basis_block[norp][pi_l][i].r[p][N[norp]/2]>0)
            {
                int sp_index=basis_block[norp][pi_l][i].r[p][N[norp]/2]-1;
                sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,0);
            }
            else
            {
                int sp_index=-basis_block[norp][pi_l][i].r[p][N[norp]/2]-1;
                sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),0);
            }
        }
        for(int q=0;q<dim_r;q++)
        {
            for(int kk=0;kk<N[norp]/2;kk++)
            {
                if(basis_block[norp][pi_r][j].r[q][kk]<0)
                {
                    sr[kk+N[norp]/2+N[norp]%2]=PQ_mat_out(P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].J,P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].pi,-P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].M,P_M[norp][-basis_block[norp][pi_r][j].r[q][kk]].PQ_Jpi_index,basis_block[norp][pi_r][j].r[q][kk]);
                }
                else
                {
                    sr[kk+N[norp]/2+N[norp]%2]=PQ_mat_out(P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].J,P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].pi,P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].M,P_M[norp][basis_block[norp][pi_r][j].r[q][kk]].PQ_Jpi_index,basis_block[norp][pi_r][j].r[q][kk]);
                }
            }
            if(N[norp]%2==1)
            {
                if(basis_block[norp][pi_r][j].r[q][N[norp]/2]>0)
                {
                    int sp_index=basis_block[norp][pi_r][j].r[q][N[norp]/2]-1;
                    sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,1);
                }
                else
                {
                    int sp_index=-basis_block[norp][pi_r][j].r[q][N[norp]/2]-1;
                    sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),1);
                }
            }
            PQ_mat_struct Q;
            if(morp==0)
            {
                Q=PQ_mat_out(2,0,2,-1,J_pm_index[norp]);
            }
            else
            {
                Q=PQ_mat_out(2,0,-2,-1,-J_pm_index[norp]);
            }
            if(N[norp]%2==0)
            {
                val[p+q*dim_l]=q_cal(sr,Q,N[norp]/2,norp);
            }
            else
            {
                val[p+q*dim_l]=q_o_cal(sr,Q,N[norp]/2,norp);
            }
            // cout<<val[p+q*dim_l]<<' ';
            PQ_mat_free(Q);
        }
        // cout<<endl;
    }
    ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp]=new double [dim_r*dim_l];
    ham_q_block[norp][pi_l][pi_r][i][j].Jpm_unnor[morp]=val;
    mat_ort_nor(ham_q_block[norp][pi_l][pi_r][i][j].Jpm[morp],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
    ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[morp]=true;
}

void ham_q_block_init()
{
    for(int norp=0;norp<2;norp++)
    {
        for(int pi_l=0;pi_l<2;pi_l++)
        {
            for(int pi_r=0;pi_r<2;pi_r++)
            {
                ham_q_block[norp][pi_l][pi_r]=new struct ham_q_block_struct *[basis_block[norp][pi_l].size()];
                for(int i=0;i<basis_block[norp][pi_l].size();i++)
                {
                    ham_q_block[norp][pi_l][pi_r][i]=new struct ham_q_block_struct [basis_block[norp][pi_r].size()];
                    for(int j=0;j<basis_block[norp][pi_r].size();j++)
                    {
                        for(int morp=0;morp<2;morp++)
                        {
                            ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp]=new double **[QQ_np_num];
                            ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[morp]=new double **[QQ_np_num];
                            ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp]=new bool *[QQ_np_num];
                            ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[morp]=false;
                            for(int k=0;k<QQ_np_num;k++)
                            {
                                ham_q_block[norp][pi_l][pi_r][i][j].Q_np[morp][k]=new double *[QQ_np_ham[k].Q_index[norp].size()];
                                ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[morp][k]=new double *[QQ_np_ham[k].Q_index[norp].size()];
                                ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp][k]=new bool [QQ_np_ham[k].Q_index[norp].size()];
                                for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                                {
                                    ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[morp][k][qi]=false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
void ham_q_block_out(int norp,int pi_l,int pi_r)
{
    
    for(int i=0;i<basis_block[norp][pi_l].size();i++)
    {
        for(int j=0;j<basis_block[norp][pi_r].size();j++)
        {
            
            for(int k=0;k<QQ_np_num;k++)
            {
                //cout<<"QQ_np k="<<k<<endl;
                for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                {
                    // if(QQ_np_ham[k].her==-1&&qi==0)
                    // {
                    //     // continue;
                    // }
                    if(basis_block[norp][pi_l][i].M==basis_block[norp][pi_r][j].M+Q_M[norp][QQ_np_ham[k].Q_index[norp][qi]-INT_MAX/2].M&&(pi_l+pi_r+Q_M[norp][QQ_np_ham[k].Q_index[norp][qi]-INT_MAX/2].pi)%2==0)
                    {
                        //cout<<"QQ_np qi="<<qi<<endl;
                        ham_q_block_insert(i,j,k,qi,0,norp,pi_l,pi_r);
                    }
                    if(basis_block[norp][pi_l][i].M==basis_block[norp][pi_r][j].M-Q_M[norp][QQ_np_ham[k].Q_index[norp][qi]-INT_MAX/2].M&&(pi_l+pi_r+Q_M[norp][QQ_np_ham[k].Q_index[norp][qi]-INT_MAX/2].pi)%2==0)
                    {
                        //cout<<"QQ_np qi="<<-qi<<endl;
                        ham_q_block_insert(i,j,k,qi,1,norp,pi_l,pi_r);
                    }
                }
            }
            if(basis_block[norp][pi_l][i].M==basis_block[norp][pi_r][j].M+2&&pi_l==pi_r)
            {
                // cout<<"J+="<<endl;
                Jpm_q_block_insert(i,j,0,norp,pi_l,pi_r);
            }
            if(basis_block[norp][pi_l][i].M==basis_block[norp][pi_r][j].M-2&&pi_l==pi_r)
            {
                // cout<<"J-="<<endl;
                Jpm_q_block_insert(i,j,1,norp,pi_l,pi_r);
            }
        }
    }
}
void ham_q_block_write()
{
    int temp_d;
    char filename[]="ham_q_block_0_0_0.bin";
    for(int norp=0;norp<2;norp++)
    {
        for(int pi_l=0;pi_l<2;pi_l++)
        {
            for(int pi_r=0;pi_r<2;pi_r++)
            {
                if(block_mat[norp][pi_l].size()<=block_store_num[norp][pi_l]&&block_mat[norp][pi_r].size()<=block_store_num[norp][pi_r])
                {
                    continue;
                }
                filename[12]='0'+norp;
                filename[14]='0'+pi_l;
                filename[16]='0'+pi_r;
                FILE *fp=fopen(filename,"wb");
                int block_num=basis_block[norp][pi_l].size()*basis_block[norp][pi_r].size();
                fwrite(&block_num,sizeof(int),1,fp);
                for(int i=0;i<basis_block[norp][pi_l].size();i++)
                {
                    for(int j=0;j<basis_block[norp][pi_r].size();j++)
                    {
                        fwrite(&(basis_block[norp][pi_l][i].M),sizeof(int),1,fp);
                        int dim_l=basis_block[norp][pi_l][i].r.size();
                        fwrite(&dim_l,sizeof(int),1,fp);
                        fwrite(&(basis_block[norp][pi_r][j].M),sizeof(int),1,fp);
                        int dim_r=basis_block[norp][pi_r][j].r.size();
                        fwrite(&dim_r,sizeof(int),1,fp);
                        for(int k=0;k<QQ_np_num;k++)
                        {
                            fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[0][k],sizeof(bool),QQ_np_ham[k].Q_index[norp].size(),fp);
                            for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                            {
                                if(ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[0][k][qi])
                                {
                                    fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Q_np[0][k][qi],sizeof(double),dim_l*dim_r,fp);
                                    fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[0][k][qi],sizeof(double),dim_l*dim_r,fp);
                                }
                            }
                            fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[1][k],sizeof(bool),QQ_np_ham[k].Q_index[norp].size(),fp);
                            for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                            {
                                if(ham_q_block[norp][pi_l][pi_r][i][j].Q_np_evalued[1][k][qi])
                                {
                                    fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Q_np[1][k][qi],sizeof(double),dim_l*dim_r,fp);
                                    fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Q_np_unnor[1][k][qi],sizeof(double),dim_l*dim_r,fp);
                                }
                            }
                            fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued,sizeof(bool),1,fp);
                            if(ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[0])
                            {
                                fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Jpm[0],sizeof(double),dim_l*dim_r,fp);
                                fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Jpm_unnor[0],sizeof(double),dim_l*dim_r,fp);
                            }
                            fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued+1,sizeof(bool),1,fp);
                            if(ham_q_block[norp][pi_l][pi_r][i][j].Jpm_evalued[1])
                            {
                                fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Jpm[1],sizeof(double),dim_l*dim_r,fp);
                                fwrite(ham_q_block[norp][pi_l][pi_r][i][j].Jpm_unnor[1],sizeof(double),dim_l*dim_r,fp);
                            }
                        }
                    }
                }
                fclose(fp);
            }
        }
    }
}
void ham_q_block_read()
{
    char filename[]="ham_q_block_0_0_0.bin";
    for(int norp=0;norp<2;norp++)
    {
        for(int pi_l=0;pi_l<2;pi_l++)
        {
            for(int pi_r=0;pi_r<2;pi_r++)
            {
                filename[12]='0'+norp;
                filename[14]='0'+pi_l;
                filename[16]='0'+pi_r;
                FILE *fp;
                if((fp=fopen(filename,"rb"))==NULL)
                {
                    continue;
                }
                int block_num;
                fread(&block_num,sizeof(int),1,fp);
                for(int p=0;p<block_num;p++)
                {
                    int M_l,dim_l,M_r,dim_r;
                    fread(&M_l,sizeof(int),1,fp);
                    fread(&dim_l,sizeof(int),1,fp);
                    fread(&M_r,sizeof(int),1,fp);
                    fread(&dim_r,sizeof(int),1,fp);
                    bool block_exsist_l=false;
                    bool block_exsist_r=false;
                    int index_l,index_r;
                    for(int i=0;i<basis_block[norp][pi_l].size();i++)
                    {
                        if(basis_block[norp][pi_l][i].M==M_l)
                        {
                            block_exsist_l=true;
                            index_l=i;
                            break;
                        }
                    }
                    for(int i=0;i<basis_block[norp][pi_r].size();i++)
                    {
                        if(basis_block[norp][pi_r][i].M==M_r)
                        {
                            block_exsist_r=true;
                            index_r=i;
                            break;
                        }
                    }
                    if(block_exsist_l&&block_exsist_r)
                    {
                        for(int k=0;k<QQ_np_num;k++)
                        {
                            fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_evalued[0][k],sizeof(bool),QQ_np_ham[k].Q_index[norp].size(),fp);
                            for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                            {
                                if(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_evalued[0][k][qi])
                                {
                                    ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np[0][k][qi]=new double [dim_l*dim_r];
                                    fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np[0][k][qi],sizeof(double),dim_l*dim_r,fp);
                                    ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_unnor[0][k][qi]=new double [dim_l*dim_r];
                                    fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_unnor[0][k][qi],sizeof(double),dim_l*dim_r,fp);
                                }
                            }

                            fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_evalued[1][k],sizeof(bool),QQ_np_ham[k].Q_index[norp].size(),fp);
                            for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                            {
                                if(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_evalued[1][k][qi])
                                {
                                    ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np[1][k][qi]=new double [dim_l*dim_r];
                                    fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np[1][k][qi],sizeof(double),dim_l*dim_r,fp);
                                    ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_unnor[1][k][qi]=new double [dim_l*dim_r];
                                    fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Q_np_unnor[1][k][qi],sizeof(double),dim_l*dim_r,fp);
                                }
                            }
                            fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_evalued,sizeof(bool),1,fp);
                            if(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_evalued[0])
                            {
                                ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm[0]=new double [dim_l*dim_r];
                                fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm[0],sizeof(double),dim_l*dim_r,fp);
                                ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_unnor[0]=new double [dim_l*dim_r];
                                fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_unnor[0],sizeof(double),dim_l*dim_r,fp);
                            }
                            fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_evalued+1,sizeof(bool),1,fp);
                            if(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_evalued[1])
                            {
                                ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm[1]=new double [dim_l*dim_r];
                                fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm[1],sizeof(double),dim_l*dim_r,fp);
                                ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_unnor[1]=new double [dim_l*dim_r];
                                fread(ham_q_block[norp][pi_l][pi_r][index_l][index_r].Jpm_unnor[1],sizeof(double),dim_l*dim_r,fp);
                            }
                        }
                    }
                    else
                    {
                        for(int k=0;k<QQ_np_num;k++)
                        {
                            bool temp_b[QQ_np_ham[k].Q_index[norp].size()];
                            fread(temp_b,sizeof(bool),QQ_np_ham[k].Q_index[norp].size(),fp);
                            for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                            {
                                if(temp_b[qi])
                                {
                                    fseek(fp,2*sizeof(double)*dim_l*dim_r,SEEK_CUR);
                                }
                            }
                            fread(temp_b,sizeof(bool),QQ_np_ham[k].Q_index[norp].size(),fp);
                            for(int qi=0;qi<QQ_np_ham[k].Q_index[norp].size();qi++)
                            {
                                if(temp_b[qi])
                                {
                                    fseek(fp,2*sizeof(double)*dim_l*dim_r,SEEK_CUR);
                                }
                            }
                            fread(temp_b,sizeof(bool),1,fp);
                            if(temp_b[0])
                            {
                                fseek(fp,2*sizeof(double)*dim_l*dim_r,SEEK_CUR);
                            }
                            fread(temp_b,sizeof(bool),1,fp);
                            if(temp_b[0])
                            {
                                fseek(fp,2*sizeof(double)*dim_l*dim_r,SEEK_CUR);
                            }
                        }
                    }
                    
                }
                fclose(fp);
            }
        }
    }
}