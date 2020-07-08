struct tran_q_block_struct
{
    double **mat;
    double **mat_unnor;
    bool *mat_evalued;
};
tran_q_block_struct **tran_q_block[2][2][2];

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

void tran_q_block_insert(int i,int j,int k,int norp,int pi_l,int pi_r)
{
    int her=Q_tran_her[norp][k];
    int Q_index=Q_M_tran_index[norp][k];
    int M_l=basis_block[norp][pi_l][i].M;
    int M_r=basis_block[norp][pi_r][j].M;
    int dim_l=basis_block[norp][pi_l][i].r.size();
    int dim_r=basis_block[norp][pi_r][j].r.size();
    int Q_k=Q_M[norp][Q_index-INT_MAX/2].J;
    double *val;
    if(tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued[k])
    {
        return;
    }
    int Rev_i=basis_block[norp][pi_l].size()-1-i;
    int Rev_j=basis_block[norp][pi_r].size()-1-j;
    if(tran_q_block[norp][pi_l][pi_r][Rev_i][Rev_j].mat_evalued[k])
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
                val[p+q*dim_l]=phase*tran_q_block[norp][pi_l][pi_r][Rev_i][Rev_j].mat_unnor[k][Rev_p+Rev_q*dim_l];
            }
        }
        tran_q_block[norp][pi_l][pi_r][i][j].mat[k]=new double [dim_r*dim_l];
        tran_q_block[norp][pi_l][pi_r][i][j].mat_unnor[k]=val;
        mat_ort_nor(tran_q_block[norp][pi_l][pi_r][i][j].mat[k],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued[k]=true;
        return;
    }
    if(tran_q_block[norp][pi_r][pi_l][j][i].mat_evalued[k])
    {
        val=new double [dim_r*dim_l];
        memcpy(val,tran_q_block[norp][pi_r][pi_l][j][i].mat_unnor[k],sizeof(double)*dim_l*dim_r);
        mkl_dimatcopy('c', 't', dim_r, dim_l, her, val, dim_r, dim_l);
        tran_q_block[norp][pi_l][pi_r][i][j].mat[k]=new double [dim_r*dim_l];
        tran_q_block[norp][pi_l][pi_r][i][j].mat_unnor[k]=val;
        mat_ort_nor(tran_q_block[norp][pi_l][pi_r][i][j].mat[k],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued[k]=true;
        return;
    }
    
    if(tran_q_block[norp][pi_r][pi_l][Rev_j][Rev_i].mat_evalued[k])
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
                val[p+q*dim_l]=phase*tran_q_block[norp][pi_r][pi_l][Rev_j][Rev_i].mat_unnor[k][Rev_q+Rev_p*dim_r];
            }
        }
        tran_q_block[norp][pi_l][pi_r][i][j].mat[k]=new double [dim_r*dim_l];
        tran_q_block[norp][pi_l][pi_r][i][j].mat_unnor[k]=val;
        mat_ort_nor(tran_q_block[norp][pi_l][pi_r][i][j].mat[k],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
        tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued[k]=true;
        return;
    }
    cout<<"block "<<i<<" in "<<basis_block[norp][pi_l].size()<<" of pi_l= "<<pi_l<<" ; operator Q_i="<<k<<"; "<<j<<" in "<<basis_block[norp][pi_r].size()<<" of pi_r= "<<pi_r<<" "<<norp<<"th neuclon"<<endl;
    PQ_mat_struct sr[(N[norp]/2+N[norp]%2)*2];
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
            Q=PQ_mat_out(Q_M[norp][Q_index-INT_MAX/2].J,Q_M[norp][Q_index-INT_MAX/2].pi,Q_M[norp][Q_index-INT_MAX/2].M,Q_M[norp][Q_index-INT_MAX/2].PQ_Jpi_index,Q_index);
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
    tran_q_block[norp][pi_l][pi_r][i][j].mat[k]=new double [dim_r*dim_l];
    tran_q_block[norp][pi_l][pi_r][i][j].mat_unnor[k]=val;
    mat_ort_nor(tran_q_block[norp][pi_l][pi_r][i][j].mat[k],val,block_mat[norp][pi_l][i].vec,block_mat[norp][pi_r][j].vec,block_mat[norp][pi_l][i].dim,block_mat[norp][pi_r][j].dim);
    tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued[k]=true;
}

void tran_q_block_init()
{
    for(int norp=0;norp<2;norp++)
    {
        for(int pi_l=0;pi_l<2;pi_l++)
        {
            for(int pi_r=0;pi_r<2;pi_r++)
            {
                tran_q_block[norp][pi_l][pi_r]=new struct tran_q_block_struct *[basis_block[norp][pi_l].size()];
                for(int i=0;i<basis_block[norp][pi_l].size();i++)
                {
                    tran_q_block[norp][pi_l][pi_r][i]=new struct tran_q_block_struct [basis_block[norp][pi_r].size()];
                    for(int j=0;j<basis_block[norp][pi_r].size();j++)
                    {
                        tran_q_block[norp][pi_l][pi_r][i][j].mat=new double *[Q_tran_num[norp]];
                        tran_q_block[norp][pi_l][pi_r][i][j].mat_unnor=new double *[Q_tran_num[norp]];
                        tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued=new bool [Q_tran_num[norp]];
                        for(int k=0;k<Q_tran_num[norp];k++)
                        {
                            tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued[k]=false;
                        }
                    }
                }
            }
        }
    }
}
void tran_q_block_out(int norp,int pi_l,int pi_r)
{
    for(int i=0;i<basis_block[norp][pi_l].size();i++)
    {
        for(int j=0;j<basis_block[norp][pi_r].size();j++)
        {
            
            for(int k=0;k<Q_tran_num[norp];k++)
            {
                if(basis_block[norp][pi_l][i].M==basis_block[norp][pi_r][j].M&&(pi_l+pi_r+Q_M[norp][Q_M_tran_index[norp][k]-INT_MAX/2].pi)%2==0)
                {
                    tran_q_block_insert(i,j,k,norp,pi_l,pi_r);
                }
            }
        }
    }
}
void tran_q_block_write()
{
    int temp_d;
    char filename[]="tran_q_block_0_0_0.bin";
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
                        fwrite(tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued,sizeof(bool),Q_tran_num[norp],fp);
                        for(int k=0;k<Q_tran_num[norp];k++)
                        {
                            if(tran_q_block[norp][pi_l][pi_r][i][j].mat_evalued[k])
                            {
                                fwrite(tran_q_block[norp][pi_l][pi_r][i][j].mat[k],sizeof(double),dim_l*dim_r,fp);
                                fwrite(tran_q_block[norp][pi_l][pi_r][i][j].mat_unnor[k],sizeof(double),dim_l*dim_r,fp);
                            }
                        }
                    }
                }
                fclose(fp);
            }
        }
    }
}
void tran_q_block_read()
{
    char filename[]="tran_q_block_0_0_0.bin";
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
                        fread(tran_q_block[norp][pi_l][pi_r][index_l][index_r].mat_evalued,sizeof(bool),Q_tran_num[norp],fp);
                        for(int k=0;k<Q_tran_num[norp];k++)
                        {
                            if(tran_q_block[norp][pi_l][pi_r][index_l][index_r].mat_evalued[k])
                            {
                                tran_q_block[norp][pi_l][pi_r][index_l][index_r].mat[k]=new double [dim_l*dim_r];
                                fread(tran_q_block[norp][pi_l][pi_r][index_l][index_r].mat[k],sizeof(double),dim_l*dim_r,fp);
                                tran_q_block[norp][pi_l][pi_r][index_l][index_r].mat_unnor[k]=new double [dim_l*dim_r];
                                fread(tran_q_block[norp][pi_l][pi_r][index_l][index_r].mat_unnor[k],sizeof(double),dim_l*dim_r,fp);
                            }
                        }
                    }
                    else
                    {
                        bool temp_b[Q_tran_num[norp]];
                        fread(temp_b,sizeof(bool),Q_tran_num[norp],fp);
                        for(int k=0;k<Q_tran_num[norp];k++)
                        {
                            if(temp_b[k])
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