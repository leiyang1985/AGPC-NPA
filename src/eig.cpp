extern "C" void dsaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv,
int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);
extern "C" void dseupd_(bool *rvec, char *howmny, int *select, double *d, double *z, int *ldz, double *sigma, char *bmat, int *n,
char *which, int *nev, double *tol, double *resid, int *ncv, double *v, int *ldv, int *iparam, int *ipntr, double *workd, double *workl,
int *lworkl, int *info);

void spe_mv(double *vout,double *vin,int pi,int norp,int k)
//vout+=spe*vin
{
    for(auto i=basis_tot_block[pi].begin();i!=basis_tot_block[pi].end();i++)
    {
        int pi_i=i->pi[0];
        int pi_j=i->pi[1];
        int ii=i->block[0];
        int jj=i->block[1];
        int dim_i=basis_block[0][pi_i][ii].r.size();
        int dim_j=basis_block[1][pi_j][jj].r.size();
        int begin=i->begin;
        int dim=i->dim;
        int end=i->end;
        for(auto p=basis_tot[pi].begin()+begin;p!=basis_tot[pi].begin()+end;p++)
        {
            int p_l=(*p)[2];
            int q_l=(*p)[5];
            for(auto q=p;q!=basis_tot[pi].begin()+end;q++)
            {
                int p_r=(*q)[2];
                int q_r=(*q)[5];
                bool p_swap=false;
                bool q_swap=false;
                if(p_l>p_r)
                {
                    swap(&p_l,&p_r);
                    p_swap=true;
                }
                if(q_l>q_r)
                {
                    swap(&q_l,&q_r);
                    q_swap=true;
                }
                if(norp==0)
                {
                    if(q_l==q_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[0][pi_i][ii].spe[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),q)]; 
                        // if(distance(basis_tot[pi].begin(),p)==0)
                        // {
                        //     cout<<pi_i<<' '<<ii<<' '<<k<<' '<<' '<<dim_i<<' '<<p_l<<' '<<p_r<<' '<<block_mat[0][pi_i][ii].spe[k][pack_k(p_l,p_r)]<<' '<<vin[distance(basis_tot[pi].begin(),q)]<<' '<<vout[distance(basis_tot[pi].begin(),p)]<<endl;
                        // }                   
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[0][pi_i][ii].spe[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),p)];
                            // if(distance(basis_tot[pi].begin(),q)==0)
                            // {
                            //     cout<<block_mat[0][pi_i][ii].spe[k][pack_k(p_l,p_r)]<<' '<<vin[distance(basis_tot[pi].begin(),p)]<<' '<<vout[distance(basis_tot[pi].begin(),q)]<<endl;
                            // }
                        }
                    }
                    
                }
                else
                {
                    if(p_l==p_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[1][pi_j][jj].spe[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        // cout<<"ham="<<norp<<' '<<distance(basis_tot[pi].begin(),p)<<' '<<distance(basis_tot[pi].begin(),q)<<' '<<jj<<' '<<block_mat[1][pi_j][jj].spe[k][pack_k(q_l,q_r)]<<endl;
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[1][pi_j][jj].spe[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                
                if(p_swap)
                {
                    swap(&p_l,&p_r);
                }
                if(q_swap)
                {
                    swap(&q_l,&q_r);
                }
            }
        }
    }
}


void PP_mv(double *vout,double *vin,int pi,int norp,int k)
//vout+=spe*vin
{
    for(auto i=basis_tot_block[pi].begin();i!=basis_tot_block[pi].end();i++)
    {
        int pi_i=i->pi[0];
        int pi_j=i->pi[1];
        int ii=i->block[0];
        int jj=i->block[1];
        int dim_i=basis_block[0][pi_i][ii].r.size();
        int dim_j=basis_block[1][pi_j][jj].r.size();
        int begin=i->begin;
        int dim=i->dim;
        int end=i->end;
        for(auto p=basis_tot[pi].begin()+begin;p!=basis_tot[pi].begin()+end;p++)
        {
            int p_l=(*p)[2];
            int q_l=(*p)[5];
            for(auto q=p;q!=basis_tot[pi].begin()+end;q++)
            {
                int p_r=(*q)[2];
                int q_r=(*q)[5];
                bool p_swap=false;
                bool q_swap=false;
                if(p_l>p_r)
                {
                    swap(&p_l,&p_r);
                    p_swap=true;
                }
                if(q_l>q_r)
                {
                    swap(&q_l,&q_r);
                    q_swap=true;
                }
                if(norp==0)
                {
                    if(q_l==q_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_i][ii].PP[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_i][ii].PP[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                else
                {
                    if(p_l==p_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_j][jj].PP[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_j][jj].PP[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                
                if(p_swap)
                {
                    swap(&p_l,&p_r);
                }
                if(q_swap)
                {
                    swap(&q_l,&q_r);
                }
            }
        }
    }
}

void PP_in_basis_mv(double *vout,double *vin,int pi,int norp,int k)
//vout+=spe*vin
{
    for(auto i=basis_tot_block[pi].begin();i!=basis_tot_block[pi].end();i++)
    {
        int pi_i=i->pi[0];
        int pi_j=i->pi[1];
        int ii=i->block[0];
        int jj=i->block[1];
        int dim_i=basis_block[0][pi_i][ii].r.size();
        int dim_j=basis_block[1][pi_j][jj].r.size();
        int begin=i->begin;
        int dim=i->dim;
        int end=i->end;
        for(auto p=basis_tot[pi].begin()+begin;p!=basis_tot[pi].begin()+end;p++)
        {
            int p_l=(*p)[2];
            int q_l=(*p)[5];
            for(auto q=p;q!=basis_tot[pi].begin()+end;q++)
            {
                int p_r=(*q)[2];
                int q_r=(*q)[5];
                bool p_swap=false;
                bool q_swap=false;
                if(p_l>p_r)
                {
                    swap(&p_l,&p_r);
                    p_swap=true;
                }
                if(q_l>q_r)
                {
                    swap(&q_l,&q_r);
                    q_swap=true;
                }
                if(norp==0)
                {
                    if(q_l==q_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_i][ii].PP_in_basis[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_i][ii].PP_in_basis[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                else
                {
                    if(p_l==p_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_j][jj].PP_in_basis[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_j][jj].PP_in_basis[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                
                if(p_swap)
                {
                    swap(&p_l,&p_r);
                }
                if(q_swap)
                {
                    swap(&q_l,&q_r);
                }
            }
        }
    }
}

void C_in_basis_mv(double *vout,double *vin,int pi,int norp,int k)
//vout+=spe*vin
{
    for(auto i=basis_tot_block[pi].begin();i!=basis_tot_block[pi].end();i++)
    {
        int pi_i=i->pi[0];
        int pi_j=i->pi[1];
        int ii=i->block[0];
        int jj=i->block[1];
        int dim_i=basis_block[0][pi_i][ii].r.size();
        int dim_j=basis_block[1][pi_j][jj].r.size();
        int begin=i->begin;
        int dim=i->dim;
        int end=i->end;
        for(auto p=basis_tot[pi].begin()+begin;p!=basis_tot[pi].begin()+end;p++)
        {
            int p_l=(*p)[2];
            int q_l=(*p)[5];
            for(auto q=p;q!=basis_tot[pi].begin()+end;q++)
            {
                int p_r=(*q)[2];
                int q_r=(*q)[5];
                bool p_swap=false;
                bool q_swap=false;
                if(p_l>p_r)
                {
                    swap(&p_l,&p_r);
                    p_swap=true;
                }
                if(q_l>q_r)
                {
                    swap(&q_l,&q_r);
                    q_swap=true;
                }
                if(norp==0)
                {
                    if(q_l==q_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_i][ii].C_in_basis[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_i][ii].C_in_basis[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                else
                {
                    if(p_l==p_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_j][jj].C_in_basis[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_j][jj].C_in_basis[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                
                if(p_swap)
                {
                    swap(&p_l,&p_r);
                }
                if(q_swap)
                {
                    swap(&q_l,&q_r);
                }
            }
        }
    }
}


void QQ_mv(double *vout,double *vin,int pi,int norp,int k)
//vout+=spe*vin
{
    for(auto i=basis_tot_block[pi].begin();i!=basis_tot_block[pi].end();i++)
    {
        int pi_i=i->pi[0];
        int pi_j=i->pi[1];
        int ii=i->block[0];
        int jj=i->block[1];
        int dim_i=basis_block[0][pi_i][ii].r.size();
        int dim_j=basis_block[1][pi_j][jj].r.size();
        int begin=i->begin;
        int dim=i->dim;
        int end=i->end;
        for(auto p=basis_tot[pi].begin()+begin;p!=basis_tot[pi].begin()+end;p++)
        {
            int p_l=(*p)[2];
            int q_l=(*p)[5];
            for(auto q=p;q!=basis_tot[pi].begin()+end;q++)
            {
                int p_r=(*q)[2];
                int q_r=(*q)[5];
                bool p_swap=false;
                bool q_swap=false;
                if(p_l>p_r)
                {
                    swap(&p_l,&p_r);
                    p_swap=true;
                }
                if(q_l>q_r)
                {
                    swap(&q_l,&q_r);
                    q_swap=true;
                }
                if(norp==0)
                {
                    if(q_l==q_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_i][ii].QQ[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_i][ii].QQ[k][pack_k(p_l,p_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                else
                {
                    if(p_l==p_r)
                    {
                        vout[distance(basis_tot[pi].begin(),p)]+=block_mat[norp][pi_j][jj].QQ[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=block_mat[norp][pi_j][jj].QQ[k][pack_k(q_l,q_r)]*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
                
                if(p_swap)
                {
                    swap(&p_l,&p_r);
                }
                if(q_swap)
                {
                    swap(&q_l,&q_r);
                }
            }
        }
    }
}

void JpJm_mv(double *vout,double *vin,int pi)
//vout+=spe*vin
{
    // cout<<"go into JpJm once"<<endl;
    double delta;
    for(auto i_l=basis_tot_block[pi].begin();i_l!=basis_tot_block[pi].end();i_l++)
    {
        int pi_i_l=i_l->pi[0];
        int pi_j_l=i_l->pi[1];
        int ii_l=i_l->block[0];
        int jj_l=i_l->block[1];
        int dim_i_l=basis_block[0][pi_i_l][ii_l].r.size();
        int dim_j_l=basis_block[1][pi_j_l][jj_l].r.size();
        int M_i_l=basis_block[0][pi_i_l][ii_l].M;
        int M_j_l=basis_block[1][pi_j_l][jj_l].M;
        int begin_l=i_l->begin;
        int dim_l=i_l->dim;
        int end_l=i_l->end;
        for(auto i_r=i_l;i_r!=basis_tot_block[pi].end();i_r++)
        {
            int pi_i_r=i_r->pi[0];
            int pi_j_r=i_r->pi[1];
            int ii_r=i_r->block[0];
            int jj_r=i_r->block[1];
            int dim_i_r=basis_block[0][pi_i_r][ii_r].r.size();
            int dim_j_r=basis_block[1][pi_j_r][jj_r].r.size();
            int M_i_r=basis_block[0][pi_i_r][ii_r].M;
            int M_j_r=basis_block[1][pi_j_r][jj_r].M;
            int begin_r=i_r->begin;
            int dim_r=i_r->dim;
            int end_r=i_r->end;
            
            if(i_l==i_r)
            {
                for(auto p=basis_tot[pi].begin()+begin_l;p!=basis_tot[pi].begin()+end_l;p++)
                {
                    int p_l=(*p)[2];
                    int q_l=(*p)[5];
                    for(auto q=p;q!=basis_tot[pi].begin()+end_r;q++)
                    {
                        int p_r=(*q)[2];
                        int q_r=(*q)[5];
                        bool p_swap=false;
                        bool q_swap=false;
                        if(p_l>p_r)
                        {
                            swap(&p_l,&p_r);
                            p_swap=true;
                        }
                        if(q_l>q_r)
                        {
                            swap(&q_l,&q_r);
                            q_swap=true;
                        }

                        if(q_l==q_r)
                        {
                            delta=block_mat[0][pi_i_l][ii_l].JpJm[pack_k(p_l,p_r)];
                            if(p_l==p_r)
                            {
                                delta+=block_mat[1][pi_j_l][jj_l].JpJm[pack_k(q_l,q_r)];
                            }
                            vout[distance(basis_tot[pi].begin(),p)]+=delta*vin[distance(basis_tot[pi].begin(),q)];
                            if(p!=q)
                            {
                                vout[distance(basis_tot[pi].begin(),q)]+=delta*vin[distance(basis_tot[pi].begin(),p)];
                            }
                        }
                        else
                        {
                            if(p_l==p_r)
                            {
                                delta=block_mat[1][pi_j_l][jj_l].JpJm[pack_k(q_l,q_r)];
                                vout[distance(basis_tot[pi].begin(),p)]+=delta*vin[distance(basis_tot[pi].begin(),q)];
                                if(p!=q)
                                {
                                    vout[distance(basis_tot[pi].begin(),q)]+=delta*vin[distance(basis_tot[pi].begin(),p)];
                                }
                            }
                        }
                        
                        if(p_swap)
                        {
                            swap(&p_l,&p_r);
                        }
                        if(q_swap)
                        {
                            swap(&q_l,&q_r);
                        }
                    }
                }
            }
            // cout<<M_i_l<<' '<<M_j_l<<' '<<M_i_r<<' '<<M_j_r<<endl;
            if(M_i_l==M_i_r+2&&pi_i_l==pi_i_r)
            {
                for(auto p=basis_tot[pi].begin()+begin_l;p!=basis_tot[pi].begin()+end_l;p++)
                {
                    int p_l=(*p)[2];
                    int q_l=(*p)[5];
                    for(auto q=basis_tot[pi].begin()+begin_r;q!=basis_tot[pi].begin()+end_r;q++)
                    {
                        int p_r=(*q)[2];
                        int q_r=(*q)[5];
                        delta=ham_q_block[0][pi_i_l][pi_i_r][ii_l][ii_r].Jpm[0][p_l+p_r*dim_i_l]*ham_q_block[1][pi_j_l][pi_j_r][jj_l][jj_r].Jpm[1][q_l+q_r*dim_j_l];
                        vout[distance(basis_tot[pi].begin(),p)]+=delta*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=delta*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
            }
            if(M_i_l==M_i_r-2&&pi_i_l==pi_i_r)
            {
                for(auto p=basis_tot[pi].begin()+begin_l;p!=basis_tot[pi].begin()+end_l;p++)
                {
                    int p_l=(*p)[2];
                    int q_l=(*p)[5];
                    for(auto q=basis_tot[pi].begin()+begin_r;q!=basis_tot[pi].begin()+end_r;q++)
                    {
                        int p_r=(*q)[2];
                        int q_r=(*q)[5];
                        delta=ham_q_block[0][pi_i_l][pi_i_r][ii_l][ii_r].Jpm[1][p_l+p_r*dim_i_l]*ham_q_block[1][pi_j_l][pi_j_r][jj_l][jj_r].Jpm[0][q_l+q_r*dim_j_l];
                        // cout<<ii_l<<' '<<ii_r<<' '<<p_l<<' '<<p_r<<' '<<ham_q_block[0][pi_i_l][pi_i_r][ii_l][ii_r].Jpm[1][p_l+p_r*dim_i_l]<<' '<<jj_l<<' '<<q_l<<' '<<q_r<<' '<<ham_q_block[1][pi_j_l][pi_j_r][jj_l][jj_r].Jpm[0][q_l+q_r*dim_j_l]<<endl;
                        vout[distance(basis_tot[pi].begin(),p)]+=delta*vin[distance(basis_tot[pi].begin(),q)];
                        if(p!=q)
                        {
                            vout[distance(basis_tot[pi].begin(),q)]+=delta*vin[distance(basis_tot[pi].begin(),p)];
                        }
                    }
                }
            }
        }
    }
}
void QQ_np_mv(double *vout,double *vin,int pi,int k)
//vout+=spe*vin
{
    int her=QQ_np_ham[k].her;
    int pi_QQ=QQ_np_ham[k].pi;
    double delta;
    for(auto i_l=basis_tot_block[pi].begin();i_l!=basis_tot_block[pi].end();i_l++)
    {
        int pi_i_l=i_l->pi[0];
        int pi_j_l=i_l->pi[1];
        int ii_l=i_l->block[0];
        int jj_l=i_l->block[1];
        int dim_i_l=basis_block[0][pi_i_l][ii_l].r.size();
        int dim_j_l=basis_block[1][pi_j_l][jj_l].r.size();
        int M_i_l=basis_block[0][pi_i_l][ii_l].M;
        int M_j_l=basis_block[1][pi_j_l][jj_l].M;
        int begin_l=i_l->begin;
        int dim_l=i_l->dim;
        int end_l=i_l->end;
        for(auto i_r=i_l;i_r!=basis_tot_block[pi].end();i_r++)
        {
            int pi_i_r=i_r->pi[0];
            int pi_j_r=i_r->pi[1];
            int ii_r=i_r->block[0];
            int jj_r=i_r->block[1];
            int dim_i_r=basis_block[0][pi_i_r][ii_r].r.size();
            int dim_j_r=basis_block[1][pi_j_r][jj_r].r.size();
            int M_i_r=basis_block[0][pi_i_r][ii_r].M;
            int M_j_r=basis_block[1][pi_j_r][jj_r].M;
            int begin_r=i_r->begin;
            int dim_r=i_r->dim;
            int end_r=i_r->end;
            for(int qi=0;qi<QQ_np_ham[k].Q_index[0].size();qi++)
            {
                int M=Q_M[0][QQ_np_ham[k].Q_index[0][qi]-INT_MAX/2].M;
                if(M_i_l==M_i_r+M&&(pi_i_l+pi_QQ+pi_i_r)%2==0)
                {
                    for(auto p=basis_tot[pi].begin()+begin_l;p!=basis_tot[pi].begin()+end_l;p++)
                    {
                        int p_l=(*p)[2];
                        int q_l=(*p)[5];
                        for(auto q=basis_tot[pi].begin()+begin_r;q!=basis_tot[pi].begin()+end_r;q++)
                        {
                            if(i_l==i_r&&distance(p,q)<0)
                            {
                                continue;
                            }
                            int p_r=(*q)[2];
                            int q_r=(*q)[5];
                            delta=her*ham_q_block[0][pi_i_l][pi_i_r][ii_l][ii_r].Q_np[0][k][qi][p_l+p_r*dim_i_l]*ham_q_block[1][pi_j_l][pi_j_r][jj_l][jj_r].Q_np[1][k][qi][q_l+q_r*dim_j_l];
                            // cout<<distance(basis_tot[pi].begin(),p)<<' '<<distance(basis_tot[pi].begin(),q)<<' '<<ham_q_block[0][pi_i_l][pi_i_r][ii_l][ii_r].Q_np[0][k][qi][p_l+p_r*dim_i_l]<<' '<<ham_q_block[1][pi_j_l][pi_j_r][jj_l][jj_r].Q_np[1][k][qi][q_l+q_r*dim_j_l]<<' '<<delta<<endl;
                            vout[distance(basis_tot[pi].begin(),p)]+=delta*vin[distance(basis_tot[pi].begin(),q)];
                            if(p!=q)
                            {
                                vout[distance(basis_tot[pi].begin(),q)]+=delta*vin[distance(basis_tot[pi].begin(),p)];
                            }
                        }
                    }
                }
                if(M!=0&&M_i_l==M_i_r-M&&(pi_i_l+pi_QQ+pi_i_r)%2==0)
                {
                    for(auto p=basis_tot[pi].begin()+begin_l;p!=basis_tot[pi].begin()+end_l;p++)
                    {
                        int p_l=(*p)[2];
                        int q_l=(*p)[5];
                        for(auto q=basis_tot[pi].begin()+begin_r;q!=basis_tot[pi].begin()+end_r;q++)
                        {
                            int p_r=(*q)[2];
                            int q_r=(*q)[5];
                            delta=her*ham_q_block[0][pi_i_l][pi_i_r][ii_l][ii_r].Q_np[1][k][qi][p_l+p_r*dim_i_l]*ham_q_block[1][pi_j_l][pi_j_r][jj_l][jj_r].Q_np[0][k][qi][q_l+q_r*dim_j_l];
                            // cout<<distance(basis_tot[pi].begin(),p)<<' '<<distance(basis_tot[pi].begin(),q)<<' '<<ham_q_block[0][pi_i_l][pi_i_r][ii_l][ii_r].Q_np[1][k][qi][p_l+p_r*dim_i_l]<<' '<<ham_q_block[1][pi_j_l][pi_j_r][jj_l][jj_r].Q_np[0][k][qi][q_l+q_r*dim_j_l]<<' '<<delta<<endl;
                            vout[distance(basis_tot[pi].begin(),p)]+=delta*vin[distance(basis_tot[pi].begin(),q)];
                            if(p!=q)
                            {
                                vout[distance(basis_tot[pi].begin(),q)]+=delta*vin[distance(basis_tot[pi].begin(),p)];
                            }
                        }
                    }
                }
            }
        }
    }
}

void h_mv_old(double *vout,double *vin,int pi)
{
    double temp[dim_tot[pi]];
    int inc=1;
    memset(vout,0,sizeof(double)*dim_tot[pi]);
    for(int norp=0;norp<2;norp++)
    {
        for(int k=0;k<sp_nlj[norp].size();k++)
        {
            memset(temp,0,sizeof(double)*dim_tot[pi]);
            spe_mv(temp,vin,pi,norp,k);
            daxpy_(dim_tot+pi,&(ob[norp][k]),temp,&inc,vout,&inc);
        }
        for(int k=0;k<PP_num[norp];k++)
        {
            memset(temp,0,sizeof(double)*dim_tot[pi]);
            PP_mv(temp,vin,pi,norp,k);
            daxpy_(dim_tot+pi,&(PP_ham[norp][k].par),temp,&inc,vout,&inc);
        }
        for(int k=0;k<QQ_num[norp];k++)
        {
            memset(temp,0,sizeof(double)*dim_tot[pi]);
            QQ_mv(temp,vin,pi,norp,k);
            daxpy_(dim_tot+pi,&(QQ_ham[norp][k].par),temp,&inc,vout,&inc);
        }
    }
    for(int k=0;k<QQ_np_num;k++)
    {
        memset(temp,0,sizeof(double)*dim_tot[pi]);
        QQ_np_mv(temp,vin,pi,k);
        daxpy_(dim_tot+pi,&(QQ_np_ham[k].par),temp,&inc,vout,&inc);
    }
}

void h_mv(double *vout,double *vin,int pi)
{
    int inc=1;
    char uplo='u';
    double alpha=1;
    double beta=1;
    char not_tran='n';
    char t_tran='t';
    memset(vout,0,sizeof(double)*dim_tot[pi]);
    for(int i=0;i<basis_tot_block[pi].size();i++)
    {
        int dim_i=basis_tot_block[pi][i].dim;
        int begin_i=basis_tot_block[pi][i].begin;
        for(int j=i;j<basis_tot_block[pi].size();j++)
        {
            int dim_j=basis_tot_block[pi][j].dim;
            int begin_j=basis_tot_block[pi][j].begin;
            if(ham_tot[pi][pack_k(i,j)].evalued)
            {
                double *h_local=ham_tot[pi][pack_k(i,j)].h;
                if(i==j)
                {
                    dsymv_(&uplo,&dim_i,&alpha,h_local,&dim_i,vin+begin_i,&inc,&alpha,vout+begin_i,&inc);
                }
                else
                {
                    dgemv_(&not_tran,&dim_i,&dim_j,&alpha,h_local,&dim_i,vin+begin_j,&inc,&alpha,vout+begin_i,&inc);
                    dgemv_(&t_tran,&dim_i,&dim_j,&alpha,h_local,&dim_i,vin+begin_i,&inc,&alpha,vout+begin_j,&inc);
                }
                // for(int ii=0;ii<dim_i;ii++)
                // {
                //     for(int jj=0;jj<dim_j;jj++)
                //     {
                //         if(isnan(h_local[ii+jj*dim_i])==1)
                //         {
                //             cout<<i<<' '<<j<<' '<<ii<<' '<<jj<<endl;
                //         }
                //     }
                // }
                // if(i==10||j==10)
                // {
                //     cin>>wat;
                // }
            }
        }
    }
}
void h_1d_out(double *h_1d,int pi)
{
    double I[dim_tot[pi]];
    memset(I,0,sizeof(double)*dim_tot[pi]);
    for(int i=0;i<dim_tot[pi];i++)
    {
        I[i]=1;
        h_mv(h_1d+i*dim_tot[pi],I,pi);
        I[i]=0;
    }
}
double JpJm_expectation(double *vin,int pi)
{
    double temp[dim_tot[pi]];
    memset(temp,0,sizeof(double)*dim_tot[pi]);
    JpJm_mv(temp,vin,pi);
    int inc=1;
    return ddot_(dim_tot+pi,vin,&inc,temp,&inc);
}
double spe_expectation(double *vin,int pi,int norp,int k)
{
    double temp[dim_tot[pi]];
    memset(temp,0,sizeof(double)*dim_tot[pi]);
    spe_mv(temp,vin,pi,norp,k);
    int inc=1;
    return ddot_(dim_tot+pi,vin,&inc,temp,&inc);
}
double PP_expectation(double *vin,int pi,int norp,int k)
{
    double temp[dim_tot[pi]];
    memset(temp,0,sizeof(double)*dim_tot[pi]);
    PP_mv(temp,vin,pi,norp,k);
    int inc=1;
    return ddot_(dim_tot+pi,vin,&inc,temp,&inc);
}
double QQ_expectation(double *vin,int pi,int norp,int k)
{
    double temp[dim_tot[pi]];
    memset(temp,0,sizeof(double)*dim_tot[pi]);
    QQ_mv(temp,vin,pi,norp,k);
    int inc=1;
    return ddot_(dim_tot+pi,vin,&inc,temp,&inc);
}
double C_in_basis_expectation(double *vin,int pi,int norp,int k)
{
    double temp[dim_tot[pi]];
    memset(temp,0,sizeof(double)*dim_tot[pi]);
    C_in_basis_mv(temp,vin,pi,norp,k);
    int inc=1;
    return ddot_(dim_tot+pi,vin,&inc,temp,&inc);
}
double PP_in_basis_expectation(double *vin,int pi,int norp,int k)
{
    double temp[dim_tot[pi]];
    memset(temp,0,sizeof(double)*dim_tot[pi]);
    PP_in_basis_mv(temp,vin,pi,norp,k);
    int inc=1;
    return ddot_(dim_tot+pi,vin,&inc,temp,&inc);
}
double QQ_np_expectation(double *vin,int pi,int k)
{
    double temp[dim_tot[pi]];
    memset(temp,0,sizeof(double)*dim_tot[pi]);
    QQ_np_mv(temp,vin,pi,k);
    int inc=1;
    return ddot_(dim_tot+pi,vin,&inc,temp,&inc);
}
int int_sta_num[2];
int *int_sta_J[2];
int *int_sta_pi[2];
int *int_sta_max[2];

void eig_input_read()
{
    FILE *fp=fopen("eig_input.dat","r");
    for(int norp=0;norp<2;norp++)
    {
        fscanf(fp,"%d",int_sta_num+norp);
        int_sta_J[norp]=new int [int_sta_num[norp]];
        int_sta_max[norp]=new int [int_sta_num[norp]];
        for(int k=0;k<int_sta_num[norp];k++)
        {
            fscanf(fp,"%d",int_sta_J[norp]+k);
        }
        for(int k=0;k<int_sta_num[norp];k++)
        {
            fscanf(fp,"%d",int_sta_max[norp]+k);
        }
    }
    fclose(fp);
}
double lapack_dsygvx( int N, int IL, int IU, double *A, double *B, double * W, double * V )
{
        int ITYPE = 1; // A*x = (lambda)*B*x
        char JOBZ = 'V'; // compute eigenvalues and eigenvectors
        char RANGE = 'I'; // compute the ILth-IUth eigenvalues
        char UPLO = 'U'; // upper triangle of A,B are stored
        //double *A = new double [N*N];
        int LDA = N; // leading dimension of A
        //double *B = new double [N*N];// will be overwritten
        int LDB = N; // leading dimension of B
        double VL, VU; // VL and VU are not referenced if RANGE='I'
        double ABSTOL = 1E-6; //absolute error tolerance of eigenvalue
        int M;// M = IU-IL+1, if RANGE = 'I'
        //double * W = new double [N];// 1st M values are eigenvalues
        double * Z = new double [N*N];// eigen vectors, the ith columncorresponds to the ith eigenvalue
        int LDZ = N;// leading dimension of Z
        int LWORK = 8*N; // length of WORK
        double * WORK = new double [LWORK];
        int * IWORK = new int [5*N];
        int * IFAIL = new int [N];
        int INFO;
        dsygvx_( &ITYPE, &JOBZ, &RANGE, &UPLO, &N, A, &LDA, B, &LDB, &VL, &VU, &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO);
        for(int i=0;i<N;i++)
        {
            cout<<W[i]<<' ';
        }
        cout<<endl;
}
void dig()
{
    FILE *fp_vec=fopen("vec.bin","wb");
    fwrite(&M_tot,sizeof(int),1,fp_vec);
    FILE *fp=fopen("eig.dat","w");
    fprintf(fp,"N[0]=%d N[1]=%d\n",N[0],N[1]);
    for(int pi=0;pi<2;pi++)
    {
        int dim=dim_tot[pi];
        int state_num=400;
        // cout<<"this is for pi="<<pi<<". how many low-lying states are you requairing:";
        // cin>>state_num;
        if(state_num>dim)
        {
            state_num=dim;
        }
        if(state_num<=0||dim<=0)
        {
            state_num=0;
            fwrite(&state_num,sizeof(int),1,fp_vec);
            fwrite(&state_num,sizeof(int),1,fp_vec);
            fflush(fp_vec);
            continue;
        }
        if((state_num<dim_tot[pi]/2&&dim_tot[pi]>300))
        // if(true)
        {
            int ido=0;
            char bmat='I';
            char which[2];
            which[0]='S';
            which[1]='A';
            int nev=state_num;
            if(state_num>dim/2)
            {
                cout<<"the state num excess half of dim"<<endl;
                nev=dim/2;
            }
            fwrite(&nev,sizeof(int),1,fp_vec);
            fwrite(&dim,sizeof(int),1,fp_vec);
            fflush(fp_vec);
            int ncv=nev*2;
            double tol=1e-9;
            double *v=new double [dim_tot[pi]*ncv];
            int ipar[11]={1,0,10000,0,0,0,1,0,0,0,0};
            int ipntr[11]={0,0,0,0,0,0,0,0,0,0,0};
            double *workd=new double [3*dim_tot[pi]];
            int lworkl=ncv*(ncv+8);
            double *workl=new double [lworkl];
            double *resid=new double [dim_tot[pi]];
            bool flag=true;
            
            int info=0;
            cout<<"dim="<<dim<<endl;
            int iter=0;
            clock_t begin=clock();
            clock_t end;
            while(flag)
            {
                dsaupd_
                (&ido, &bmat, &dim, which, &nev, &tol, resid, &ncv, v, &dim, ipar, ipntr, workd, workl, &lworkl, &info);
                if(iter%20==0)
                {
                    end=clock();
                    cout<<"arpack iteration: "<<iter<<"th iteration with ido="<<ido<<" time gap="<<double(end-begin)/CLOCKS_PER_SEC<<endl;
                    begin=end;
                }
                iter++;
                if(ido==-1||ido==1)
                {
                    // JpJm_mv(workd+ipntr[1]-1,workd+ipntr[0]-1,pi);
                    h_mv(workd+ipntr[1]-1,workd+ipntr[0]-1,pi);
                }
                // else if(ido==2)
                // {
                //     ove_mv(workd+ipntr[1]-1,workd+ipntr[0]-1,pi);
                // }
                else
                {
                    flag=false;
                }
            }
            cout<<"finish iteration. Go on then:"<<endl;
            // cin>>wat;
            bool rvec=true;
            char howmny='A';
            double sigma=0;
            int *select=new int [ncv];
            double *e=new double [ncv*2];
            double *z=new double [dim*nev];
            dseupd_(&rvec,&howmny,select,e,z,&nev,&sigma, &bmat, &dim, which,&nev,&tol,resid,&ncv,v,&dim, ipar, ipntr, workd,workl,&lworkl,&info);
            
            vector<pair<int,int>> J_order;
            cout<<"arpack="<<endl;
            for(int i=0;i<nev;i++)
            {
                // cout<<"J^pi= ";
                double exp_JpJm=JpJm_expectation(v+i*dim,pi);
                int J=abs(round(sqrt(4*(exp_JpJm+M_tot*M_tot*0.25-M_tot*0.5)+1)-1));
                bool isJexist=false;
                int J_order_current;
                bool isPrint=false;
                for(int k=0;k<J_order.size();k++)
                {
                    if(J_order[k].first==J)
                    {
                        J_order[k].second++;
                        J_order_current=J_order[k].second;
                        // cout<<J<<"_"<<J_order[k].second<<'^'<<pi<<' ';
                        isJexist=true;
                        fwrite(&J,sizeof(int),1,fp_vec);
                        fwrite(&(J_order[k].second),sizeof(int),1,fp_vec);
                        fwrite(&pi,sizeof(int),1,fp_vec);
                        fwrite(e+i,sizeof(double),1,fp_vec);
                        fwrite(v+i*dim,sizeof(double),dim,fp_vec);
                        fflush(fp_vec);
                        break;
                    }
                }
                // isPrint=true;
                if(!isJexist)
                {
                    pair<int,int> temp;
                    temp.first=J;
                    temp.second=1;
                    J_order_current=1;
                    J_order.emplace_back(temp);
                    // cout<<J<<"_"<<1<<'^'<<pi<<' ';
                    fwrite(&J,sizeof(int),1,fp_vec);
                    fwrite(&(temp.second),sizeof(int),1,fp_vec);
                    fwrite(&pi,sizeof(int),1,fp_vec);
                    fwrite(e+i,sizeof(double),1,fp_vec);
                    fwrite(v+i*dim,sizeof(double),dim,fp_vec);
                    fflush(fp_vec);
                }
                for(int k=0;k<int_sta_num[pi];k++)
                {
                    if(int_sta_J[pi][k]==J&&int_sta_max[pi][k]>=J_order_current)
                    {
                        isPrint=true;
                        break;
                    }
                }
                double delta;
                if(isPrint)
                {
                    if(fabs(round(sqrt(4*(exp_JpJm+M_tot*M_tot*0.25-M_tot*0.5)+1)-1)-sqrt(4*(exp_JpJm+M_tot*M_tot*0.25-M_tot*0.5)+1)+1)>0.01)
                    {
                        cout<<"non interger J, maybe energy degenercy or error"<<endl;
                    }
                    cout<<"J^pi="<<J<<"_"<<J_order_current<<"^"<<pi<<" E= "<<e[i]<<" JpJm= "<<exp_JpJm<<' ';
                    fprintf(fp,"J^pi= %d_%d^%d E= %f <JpJm>= %f ",J,J_order_current,pi,e[i],exp_JpJm);
                    if(yorn_eig_display[0]=='y')
                    {
                        cout<<"CP_num= ";
                        fprintf(fp,"CP_num= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<sp_nlj[norp].size();k++)
                            {
                                if(N[norp]%2==1)
                                {
                                    delta=C_in_basis_expectation(v+i*dim,pi,norp,k);
                                    fprintf(fp,"%f ",delta);
                                    cout<<delta<<' ';
                                }
                            }
                            for(int k=0;k<P_Jpi_in_basis_num[norp];k++)
                            {
                                delta=PP_in_basis_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[1]=='y')
                    {
                        cout<<"sp_occ= ";
                        fprintf(fp,"sp_occ= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<sp_nlj[norp].size();k++)
                            {
                                delta=spe_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[2]=='y')
                    {
                        cout<<"<PP>= ";
                        fprintf(fp,"<PP>= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<PP_num[norp];k++)
                            {
                                delta=PP_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[3]=='y')
                    {
                        cout<<"<QQ>= ";
                        fprintf(fp,"<QQ>= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<QQ_num[norp];k++)
                            {
                                delta=QQ_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[4]=='y')
                    {
                        cout<<"<QQ_np>= ";
                         fprintf(fp,"<QQ_np>= ");
                        for(int k=0;k<QQ_np_num;k++)
                        {
                            delta=QQ_np_expectation(v+i*dim,pi,k);
                            fprintf(fp,"%f ",delta);
                            cout<<delta<<' ';
                        }
                    }
                    cout<<endl;
                    fprintf(fp,"\n");
                }
            }
        }
        else
        {
            fwrite(&state_num,sizeof(int),1,fp_vec);
            fwrite(&dim,sizeof(int),1,fp_vec);
            fflush(fp_vec);
            double *h_1d=new double [dim*dim];
            memset(h_1d,0,sizeof(double)*dim*dim);
            h_1d_out(h_1d,pi);
            // for(int i=0;i<dim;i++)
            // {
            //     for(int j=0;j<dim;j++)
            //     {
            //         if(isnan(h_1d[i+j*dim])==1)
            //         {
            //             cout<<i<<' '<<j<<' '<<h_1d[i+j*dim]<<endl;
            //             cin>>wat;
            //         }
            //     }
            //     // cout<<endl;
            // }
            char jobz='V';
            char range='I';
            char uplo='U';
            double *work=new double [1];
            int lwork=-1;
            double vl,vu;
            int il=1;
            int iu=state_num;
            int m;
            double abstol=1e-14;
            int *iwork=new int [1];
            int liwork=-1;
            double e[state_num];
            double *v=new double [state_num*dim_tot[pi]];
            int isuppz[2*state_num];
            int info;
            dsyevr_(&jobz,&range,&uplo,dim_tot+pi,h_1d,dim_tot+pi, &vl, &vu, &il, &iu, &abstol, &m, e, v, dim_tot+pi,isuppz, work,&lwork, iwork, &liwork, &info);
            lwork=work[0];
            delete [] work;
            work=new double [lwork];
            liwork=iwork[0];
            delete [] iwork;
            iwork=new int [liwork];
            dsyevr_(&jobz,&range,&uplo,dim_tot+pi,h_1d,dim_tot+pi, &vl, &vu, &il, &iu, &abstol, &m, e, v, dim_tot+pi,isuppz, work,&lwork, iwork, &liwork, &info);
            delete [] work;
            delete [] iwork;
            // for(int i=0;i<8;i++)
            // {
            //     cout<<e[i]<<' ';
            // }
            // cout<<endl;
            vector<pair<int,int>> J_order;
            cout<<"dsyevr="<<endl;
            for(int i=0;i<state_num;i++)
            {
                // cout<<"J^pi= ";
                double exp_JpJm=JpJm_expectation(v+i*dim,pi);
                int J=abs(round(sqrt(4*(exp_JpJm+M_tot*M_tot*0.25-M_tot*0.5)+1)-1));
                
                bool isJexist=false;
                int J_order_current;
                bool isPrint=false;
                for(int k=0;k<J_order.size();k++)
                {
                    if(J_order[k].first==J)
                    {
                        J_order[k].second++;
                        J_order_current=J_order[k].second;
                        // cout<<J<<"_"<<J_order[k].second<<'^'<<pi<<' ';
                        isJexist=true;
                        fwrite(&J,sizeof(int),1,fp_vec);
                        fwrite(&(J_order[k].second),sizeof(int),1,fp_vec);
                        fwrite(&pi,sizeof(int),1,fp_vec);
                        fwrite(e+i,sizeof(double),1,fp_vec);
                        fwrite(v+i*dim,sizeof(double),dim,fp_vec);
                        fflush(fp_vec);
                        break;
                    }
                }
                // isPrint=true;
                if(!isJexist)
                {
                    pair<int,int> temp;
                    temp.first=J;
                    temp.second=1;
                    J_order_current=1;
                    J_order.emplace_back(temp);
                    // cout<<J<<"_"<<1<<'^'<<pi<<' ';
                    fwrite(&J,sizeof(int),1,fp_vec);
                    fwrite(&(temp.second),sizeof(int),1,fp_vec);
                    fwrite(&pi,sizeof(int),1,fp_vec);
                    fwrite(e+i,sizeof(double),1,fp_vec);
                    fwrite(v+i*dim,sizeof(double),dim,fp_vec);
                    fflush(fp_vec);
                }
                for(int k=0;k<int_sta_num[pi];k++)
                {
                    if(int_sta_J[pi][k]==J&&int_sta_max[pi][k]>=J_order_current)
                    {
                        isPrint=true;
                        break;
                    }
                }
                // isPrint=true;
                double delta;
                if(isPrint)
                {
                    if(fabs(round(sqrt(4*(exp_JpJm+M_tot*M_tot*0.25-M_tot*0.5)+1)-1)-sqrt(4*(exp_JpJm+M_tot*M_tot*0.25-M_tot*0.5)+1)+1)>0.01)
                    {
                        cout<<"non interger J, maybe energy degenercy or error"<<endl;
                    }
                    cout<<"J^pi="<<J<<"_"<<J_order_current<<"^"<<pi<<" E= "<<e[i]<<" JpJm= "<<exp_JpJm<<' ';
                    fprintf(fp,"J^pi= %d_%d^%d E= %f <JpJm>= %f ",J,J_order_current,pi,e[i],exp_JpJm);
                    if(yorn_eig_display[0]=='y')
                    {
                        cout<<"CP_num= ";
                        fprintf(fp,"CP_num= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<sp_nlj[norp].size();k++)
                            {
                                if(N[norp]%2==1)
                                {
                                    delta=C_in_basis_expectation(v+i*dim,pi,norp,k);
                                    fprintf(fp,"%f ",delta);
                                    cout<<delta<<' ';
                                }
                            }
                            for(int k=0;k<P_Jpi_in_basis_num[norp];k++)
                            {
                                delta=PP_in_basis_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[1]=='y')
                    {
                        cout<<"sp_occ= ";
                        fprintf(fp,"sp_occ= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<sp_nlj[norp].size();k++)
                            {
                                delta=spe_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[2]=='y')
                    {
                        cout<<"<PP>= ";
                        fprintf(fp,"<PP>= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<PP_num[norp];k++)
                            {
                                delta=PP_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[3]=='y')
                    {
                        cout<<"<QQ>= ";
                        fprintf(fp,"<QQ>= ");
                        for(int norp=0;norp<2;norp++)
                        {
                            for(int k=0;k<QQ_num[norp];k++)
                            {
                                delta=QQ_expectation(v+i*dim,pi,norp,k);
                                fprintf(fp,"%f ",delta);
                                cout<<delta<<' ';
                            }
                        }
                    }
                    if(yorn_eig_display[4]=='y')
                    {
                        cout<<"<QQ_np>= ";
                         fprintf(fp,"<QQ_np>= ");
                        for(int k=0;k<QQ_np_num;k++)
                        {
                            delta=QQ_np_expectation(v+i*dim,pi,k);
                            fprintf(fp,"%f ",delta);
                            cout<<delta<<' ';
                        }
                    }
                    cout<<endl;
                    fprintf(fp,"\n");
                }
            }
        }
        
    }
    fclose(fp);
    fclose(fp_vec);
}