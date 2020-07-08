struct ham_tot_block_struct
{
    bool evalued;
    double *h;
};
struct ham_tot_block_struct *ham_tot[2];//0 for parity 0 1 for parity 1
void ham_tot_init()
//creat ham_tot[] array to store the total H matrix block with packed form
{
    for(int pi_tot=0;pi_tot<=1;pi_tot++)
    {
        int block_dim=basis_tot_block[pi_tot].size();
        ham_tot[pi_tot]=new struct ham_tot_block_struct [block_dim*(block_dim+1)/2];
        for(int i=0;i<block_dim;i++)
        {
            for(int j=i;j<block_dim;j++)
            {
                ham_tot[pi_tot][pack_k(i,j)].evalued=false;
            }
        }
    }
    
}
void spe_to_ham_tot_evalue()
//evalue ham_tot array with sp matrix, it should be called first before the PP_to_ and QQ_to_
{
    int inc=1;
    for(int pi_tot=0;pi_tot<2;pi_tot++)
    {
         int block_dim=basis_tot_block[pi_tot].size();
        for(int i=0;i<block_dim;i++)
        {
            // cout<<"i="<<i<<endl;
            ham_tot[pi_tot][pack_k(i,i)].evalued=true;
            int dim_i=basis_tot_block[pi_tot][i].dim;
            int dim2_i=dim_i*dim_i;
            ham_tot[pi_tot][pack_k(i,i)].h=new double [dim_i*dim_i];
            memset(ham_tot[pi_tot][pack_k(i,i)].h,0,sizeof(double)*dim_i*dim_i);
            double *h_local=ham_tot[pi_tot][pack_k(i,i)].h;
            int basis_tot_begin=basis_tot_block[pi_tot][i].begin;
            int basis_tot_end=basis_tot_block[pi_tot][i].end;
            double *h_temp=new double [dim_i*dim_i];
            // cout<<"dim_i="<<dim_i<<endl;
            // cout<<"begin="<<basis_tot_begin<<endl;
            // cout<<"end="<<basis_tot_end<<endl;
            for(int norp=0;norp<2;norp++)
            {
                // cout<<"norp="<<norp<<endl;
                int pi_local=basis_tot_block[pi_tot][i].pi[norp];
                int block_local=basis_tot_block[pi_tot][i].block[norp];
                for(int k=0;k<sp_nlj[norp].size();k++)
                {
                    // cout<<"k="<<k<<endl;
                    if(fabs(ob[norp][k])>1e-7)
                    {
                        memset(h_temp,0,sizeof(double)*dim_i*dim_i);
                        double *spe_local=block_mat[norp][pi_local][block_local].spe[k];
                        for(int p=basis_tot_begin;p<basis_tot_end;p++)
                        {
                            // cout<<"p="<<p<<endl;
                            int p_local;
                            int p_other;
                            if(norp==0)
                            {
                                p_local=basis_tot[pi_tot][p][2];
                                p_other=basis_tot[pi_tot][p][5];
                            }
                            else
                            {
                                p_local=basis_tot[pi_tot][p][5];
                                p_other=basis_tot[pi_tot][p][2];
                            }
                            for(int q=p;q<basis_tot_end;q++)
                            {
                                // cout<<"q="<<q<<endl;
                                int q_local;
                                if(norp==0)
                                {
                                    if(p_other!=basis_tot[pi_tot][q][5])
                                    {
                                        continue;
                                    }
                                    q_local=basis_tot[pi_tot][q][2];
                                }
                                else
                                {
                                    if(p_other!=basis_tot[pi_tot][q][2])
                                    {
                                        continue;
                                    }
                                    q_local=basis_tot[pi_tot][q][5];
                                }
                                if(p_local<=q_local)
                                {
                                    h_temp[p-basis_tot_begin+(q-basis_tot_begin)*dim_i]=spe_local[pack_k(p_local,q_local)];
                                }
                                else
                                {
                                    h_temp[p-basis_tot_begin+(q-basis_tot_begin)*dim_i]=spe_local[pack_k(q_local,p_local)];
                                }
                            }
                        }
                        daxpy(&dim2_i,&(ob[norp][k]),h_temp,&inc,h_local,&inc);
                    }
                }
            }          
            delete [] h_temp;       
        }
    }
}

void PP_to_ham_tot_evalue()
//evalue ham_tot array with PP matrix, it should be called after the spe_to_
{

    int inc=1;
    for(int pi_tot=0;pi_tot<2;pi_tot++)
    {
        int block_dim=basis_tot_block[pi_tot].size();
        for(int i=0;i<block_dim;i++)
        {
            int dim_i=basis_tot_block[pi_tot][i].dim;
            int dim2_i=dim_i*dim_i;
            double *h_local=ham_tot[pi_tot][pack_k(i,i)].h;
            int basis_tot_begin=basis_tot_block[pi_tot][i].begin;
            int basis_tot_end=basis_tot_block[pi_tot][i].end;
            double *h_temp=new double [dim_i*dim_i];
            for(int norp=0;norp<2;norp++)
            {
                int pi_local=basis_tot_block[pi_tot][i].pi[norp];
                int block_local=basis_tot_block[pi_tot][i].block[norp];
                for(int k=0;k<PP_num[norp];k++)
                {
                    double par=PP_ham[norp][k].par;
                    if(fabs(par)>1e-7)
                    {
                        memset(h_temp,0,sizeof(double)*dim_i*dim_i);
                        double *PP_local=block_mat[norp][pi_local][block_local].PP[k];
                        for(int p=basis_tot_begin;p<basis_tot_end;p++)
                        {
                            int p_local;
                            int p_other;
                            if(norp==0)
                            {
                                p_local=basis_tot[pi_tot][p][2];
                                p_other=basis_tot[pi_tot][p][5];
                            }
                            else
                            {
                                p_local=basis_tot[pi_tot][p][5];
                                p_other=basis_tot[pi_tot][p][2];
                            }
                            for(int q=p;q<basis_tot_end;q++)
                            {
                                // cout<<"q="<<q<<endl;
                                int q_local;
                                if(norp==0)
                                {
                                    if(p_other!=basis_tot[pi_tot][q][5])
                                    {
                                        continue;
                                    }
                                    q_local=basis_tot[pi_tot][q][2];
                                }
                                else
                                {
                                    if(p_other!=basis_tot[pi_tot][q][2])
                                    {
                                        continue;
                                    }
                                    q_local=basis_tot[pi_tot][q][5];
                                }
                                if(p_local<=q_local)
                                {
                                    h_temp[p-basis_tot_begin+(q-basis_tot_begin)*dim_i]=PP_local[pack_k(p_local,q_local)];
                                }
                                else
                                {
                                    h_temp[p-basis_tot_begin+(q-basis_tot_begin)*dim_i]=PP_local[pack_k(q_local,p_local)];
                                }
                            }
                        }
                        daxpy(&dim2_i,&par,h_temp,&inc,h_local,&inc);
                    }
                }
            }          
            delete [] h_temp;
        }
    }
}

void QQ_to_ham_tot_evalue()
//evalue ham_tot array with QQ matrix, it should be called after the spe_to_
{

    int inc=1;
    for(int pi_tot=0;pi_tot<2;pi_tot++)
    {
        int block_dim=basis_tot_block[pi_tot].size();
        for(int i=0;i<block_dim;i++)
        {
            int dim_i=basis_tot_block[pi_tot][i].dim;
            int dim2_i=dim_i*dim_i;
            double *h_local=ham_tot[pi_tot][pack_k(i,i)].h;
            int basis_tot_begin=basis_tot_block[pi_tot][i].begin;
            int basis_tot_end=basis_tot_block[pi_tot][i].end;
            double *h_temp=new double [dim_i*dim_i];
            for(int norp=0;norp<2;norp++)
            {
                int pi_local=basis_tot_block[pi_tot][i].pi[norp];
                int block_local=basis_tot_block[pi_tot][i].block[norp];
                for(int k=0;k<QQ_num[norp];k++)
                {
                    double par=QQ_ham[norp][k].par;
                    if(fabs(par)>1e-7)
                    {
                        memset(h_temp,0,sizeof(double)*dim_i*dim_i);
                        double *QQ_local=block_mat[norp][pi_local][block_local].QQ[k];
                        for(int p=basis_tot_begin;p<basis_tot_end;p++)
                        {
                            int p_local;
                            int p_other;
                            if(norp==0)
                            {
                                p_local=basis_tot[pi_tot][p][2];
                                p_other=basis_tot[pi_tot][p][5];
                            }
                            else
                            {
                                p_local=basis_tot[pi_tot][p][5];
                                p_other=basis_tot[pi_tot][p][2];
                            }
                            for(int q=p;q<basis_tot_end;q++)
                            {
                                // cout<<"q="<<q<<endl;
                                int q_local;
                                if(norp==0)
                                {
                                    if(p_other!=basis_tot[pi_tot][q][5])
                                    {
                                        continue;
                                    }
                                    q_local=basis_tot[pi_tot][q][2];
                                }
                                else
                                {
                                    if(p_other!=basis_tot[pi_tot][q][2])
                                    {
                                        continue;
                                    }
                                    q_local=basis_tot[pi_tot][q][5];
                                }
                                if(p_local<=q_local)
                                {
                                    h_temp[p-basis_tot_begin+(q-basis_tot_begin)*dim_i]=QQ_local[pack_k(p_local,q_local)];
                                }
                                else
                                {
                                    h_temp[p-basis_tot_begin+(q-basis_tot_begin)*dim_i]=QQ_local[pack_k(q_local,p_local)];
                                }
                            }
                        }
                        daxpy(&dim2_i,&par,h_temp,&inc,h_local,&inc);
                    }
                }
            }          
            delete [] h_temp;
        }
    }
}

void QQ_np_to_ham_tot_evalue()
//evalue ham_tot array with QQ matrix, it should be called after the spe_to_
{

    int inc=1;
    for(int pi_tot=0;pi_tot<2;pi_tot++)
    {
        int block_dim=basis_tot_block[pi_tot].size();
        for(int i=0;i<block_dim;i++)
        {
            int dim_i=basis_tot_block[pi_tot][i].dim;
            int M_i_0=basis_tot_block[pi_tot][i].M[0];
            int M_i_1=basis_tot_block[pi_tot][i].M[1];
            int pi_i_0=basis_tot_block[pi_tot][i].pi[0];
            int pi_i_1=basis_tot_block[pi_tot][i].pi[1];
            int block_i_0=basis_tot_block[pi_tot][i].block[0];
            int block_i_1=basis_tot_block[pi_tot][i].block[1];
            int dim_i_0=basis_block[0][pi_i_0][block_i_0].r.size();
            int dim_i_1=basis_block[1][pi_i_1][block_i_1].r.size();
            int i_begin=basis_tot_block[pi_tot][i].begin;
            int i_end=basis_tot_block[pi_tot][i].end;
            for(int j=i;j<block_dim;j++)
            {
                int dim_j=basis_tot_block[pi_tot][j].dim;
                int M_j_0=basis_tot_block[pi_tot][j].M[0];
                int M_j_1=basis_tot_block[pi_tot][j].M[1];
                int pi_j_0=basis_tot_block[pi_tot][j].pi[0];
                int pi_j_1=basis_tot_block[pi_tot][j].pi[1];
                int block_j_0=basis_tot_block[pi_tot][j].block[0];
                int block_j_1=basis_tot_block[pi_tot][j].block[1];
                int dim_j_0=basis_block[0][pi_j_0][block_j_0].r.size();
                int dim_j_1=basis_block[1][pi_j_1][block_j_1].r.size();
                int j_begin=basis_tot_block[pi_tot][j].begin;
                int j_end=basis_tot_block[pi_tot][j].end;
                int dim2=dim_i*dim_j;
                double *h_local;
                double *h_temp=new double [dim2];
                if(ham_tot[pi_tot][pack_k(i,j)].evalued)
                {
                    h_local=ham_tot[pi_tot][pack_k(i,j)].h;
                }
                for(int k=0;k<QQ_np_num;k++)
                {
                    double par=QQ_np_ham[k].par;
                    int her=QQ_np_ham[k].her;
                    int pi=QQ_np_ham[k].pi;
                    par*=her;
                    if(fabs(par)>1e-7)
                    {
                        for(int qi=0;qi<QQ_np_ham[k].Q_index[0].size();qi++)
                        {
                            int M=Q_M[0][QQ_np_ham[k].Q_index[0][qi]-INT_MAX/2].M;
                            if(M_i_0==M_j_0+M&&(pi_i_0+pi+pi_j_0)%2==0)
                            {
                                // if(i==10&&j==10)
                                    // {
                                    //     cin>>wat;
                                    // }
                                if(!ham_tot[pi_tot][pack_k(i,j)].evalued)
                                {
                                    // if(i==10&&j==10)
                                    // {
                                    //     cin>>wat;
                                    // }
                                    ham_tot[pi_tot][pack_k(i,j)].evalued=true;
                                    ham_tot[pi_tot][pack_k(i,j)].h=new double [dim2];
                                    h_local=ham_tot[pi_tot][pack_k(i,j)].h;
                                    memset(h_local,0,sizeof(double)*dim2);
                                }
                                memset(h_temp,0,sizeof(double)*dim2);
                                for(int p=i_begin;p<i_end;p++)
                                {
                                    int p_0=basis_tot[pi_tot][p][2];
                                    int p_1=basis_tot[pi_tot][p][5];
                                    for(auto q=j_begin;q<j_end;q++)
                                    {
                                        if(i==j&&p>q)
                                        {
                                            continue;
                                        }
                                        int q_0=basis_tot[pi_tot][q][2];
                                        int q_1=basis_tot[pi_tot][q][5];
                                        h_temp[p-i_begin+(q-j_begin)*dim_i]=ham_q_block[0][pi_i_0][pi_j_0][block_i_0][block_j_0].Q_np[0][k][qi][p_0+q_0*dim_i_0]*ham_q_block[1][pi_i_1][pi_j_1][block_i_1][block_j_1].Q_np[1][k][qi][p_1+q_1*dim_i_1];
                                        // if(isnan(h_temp[p-i_begin+(q-j_begin)*dim_i])==1)
                                        // {
                                        //     cin>>wat;
                                        // }
                                    }
                                }
                                daxpy(&dim2,&par,h_temp,&inc,h_local,&inc);
                            }
                            if(M!=0&&M_i_0==M_j_0-M&&(pi_i_0+pi+pi_j_0)%2==0)
                            {
                                if(!ham_tot[pi_tot][pack_k(i,j)].evalued)
                                {
                                    // if(i==10&&j==10)
                                    // {
                                    //     cin>>wat;
                                    // }
                                    ham_tot[pi_tot][pack_k(i,j)].evalued=true;
                                    ham_tot[pi_tot][pack_k(i,j)].h=new double [dim2];
                                    h_local=ham_tot[pi_tot][pack_k(i,j)].h;
                                    memset(h_local,0,sizeof(double)*dim2);
                                }
                                memset(h_temp,0,sizeof(double)*dim2);
                                for(int p=i_begin;p<i_end;p++)
                                {
                                    int p_0=basis_tot[pi_tot][p][2];
                                    int p_1=basis_tot[pi_tot][p][5];
                                    for(auto q=j_begin;q<j_end;q++)
                                    {
                                        int q_0=basis_tot[pi_tot][q][2];
                                        int q_1=basis_tot[pi_tot][q][5];
                                        h_temp[p-i_begin+(q-j_begin)*dim_i]=ham_q_block[0][pi_i_0][pi_j_0][block_i_0][block_j_0].Q_np[1][k][qi][p_0+q_0*dim_i_0]*ham_q_block[1][pi_i_1][pi_j_1][block_i_1][block_j_1].Q_np[0][k][qi][p_1+q_1*dim_i_1];
                                        // if(isnan(h_temp[p-i_begin+(q-j_begin)*dim_i])==1)
                                        // {
                                        //     cin>>wat;
                                        // }
                                    }
                                }
                                daxpy(&dim2,&par,h_temp,&inc,h_local,&inc);
                            }
                        }
                    }
                }          
                delete [] h_temp;
            }
        }
    }
}

