#define mm_limit 6
double *mm_out(int *ser,int &norp);
double *mm_cal(int *&ser,int &norp)
{
    int size=ser[0];
    if(size==4)
    {
        double *mat_l;
        double *mat_r;
        bool T_l=false,T_r=false;
        if(abs(ser[1])<INT_MAX/2)
        {
            if(ser[1]>=0)
            {
                mat_l=pow_gamma[norp][ser[1]][ser[2]];
            }
            else
            {
                mat_l=pow_mgamma[norp][-ser[1]][ser[2]];
            }
            if(ser[2]%2==0)
            {
                T_l=true;
            }
        }
        else
        {
            if(ser[1]>0)
            {
                mat_l=pow_q[norp][ser[1]-INT_MAX/2][ser[2]];
            }
            else
            {
                mat_l=pow_qT[norp][-ser[1]-INT_MAX/2][ser[2]];
            } 
        }
        if(abs(ser[3])<INT_MAX/2)
        {
            if(ser[3]>=0)
            {
                mat_r=pow_gamma[norp][ser[3]][ser[4]];
            }
            else
            {
                mat_r=pow_mgamma[norp][-ser[3]][ser[4]];
            }
            if(ser[4]%2==0)
            {
                T_r=true;
            }
        }
        else
        {
            if(ser[3]>0)
            {
                mat_r=pow_q[norp][ser[3]-INT_MAX/2][ser[4]];
            }
            else
            {
                mat_r=pow_qT[norp][-ser[3]-INT_MAX/2][ser[4]];
            }
            
                      
        }
        if(T_l)
        {
            char side='l';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            double * mat_out=new double [sp_nljm_dim2[norp]];
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            return mat_out;
        }
        if(T_r)
        {
            char side='r';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            double * mat_out=new double [sp_nljm_dim2[norp]];
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_r,&dim,mat_l,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            return mat_out;
        }
        char opa='n';
        char opb='n';
        int dim=sp_nljm[norp].size();
        double alpha=1;
        double beta=0;
        double * mat_out=new double [sp_nljm_dim2[norp]];
        // t3=currentcycles();
        dgemm_(&opa,&opb,&dim,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
        // t4=currentcycles();
        // t_sum+=t4-t3;
        return mat_out;
    }
    else if(size==6)
    {
        if(abs(ser[1])<INT_MAX/2&&ser[2]%2==0)
        {
            int ser_r[5];
            ser_r[0]=4;
            memcpy(ser_r+1,ser+3,sizeof(int)*(size-2));
            double *mat_r=mm_out(ser_r,norp);
            double *mat_l;
            if(ser[1]>=0)
            {
                mat_l=pow_gamma[norp][ser[1]][ser[2]];
            }
            else
            {
                mat_l=pow_mgamma[norp][-ser[1]][ser[2]];
            }
            double *mat_out=new double [sp_nljm_dim2[norp]];
            char side='l';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(ser_r[0]>=mm_limit)
            {
                delete [] mat_r;
            }
            return mat_out;
        }
        if(abs(ser[5])<INT_MAX/2&&ser[6]%2==0)
        {
            int ser_l[5];
            ser_l[0]=4;
            memcpy(ser_l+1,ser+1,sizeof(int)*4);
            double *mat_r;
            double *mat_l=mm_out(ser_l,norp);
            if(ser[5]>=0)
            {
                mat_r=pow_gamma[norp][ser[5]][ser[6]];
            }
            else
            {
                mat_r=pow_mgamma[norp][-ser[5]][ser[6]];
            }
            double *mat_out=new double [sp_nljm_dim2[norp]];
            char side='r';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_r,&dim,mat_l,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(ser_l[0]>=mm_limit)
            {
                delete [] mat_l;
            }
            return mat_out;
        }
        if(abs(ser[1])>=INT_MAX/2&&ser[1]==-ser[3]&&ser[2]==ser[4])
        {
            int ser_l[5];
            ser_l[0]=4;
            memcpy(ser_l+1,ser+1,sizeof(int)*4);
            double *mat_r;
            double *mat_l=mm_out(ser_l,norp);
            if(abs(ser[5])<INT_MAX/2)
            {
                if(ser[5]>=0)
                {
                    mat_r=pow_gamma[norp][ser[5]][ser[6]];
                }
                else
                {
                    mat_r=pow_mgamma[norp][-ser[5]][ser[6]];
                }
            }
            else
            {
                if(ser[5]>0)
                {
                    mat_r=pow_q[norp][ser[5]-INT_MAX/2][ser[6]];
                }
                else
                {
                    mat_r=pow_qT[norp][-ser[5]-INT_MAX/2][ser[6]];
                }
                
            }
            double *mat_out=new double [sp_nljm_dim2[norp]];
            char side='l';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(ser_l[0]>=mm_limit)
            {
                delete [] mat_l;
            }
            return mat_out;
        }
        if(abs(ser[3])>=INT_MAX/2&&ser[3]==-ser[5]&&ser[4]==ser[6])
        {
            int ser_r[5];
            ser_r[0]=4;
            memcpy(ser_r+1,ser+3,sizeof(int)*4);
            double *mat_r=mm_out(ser_r,norp);
            double *mat_l;
            if(abs(ser[1])<INT_MAX/2)
            {
                if(ser[1]>=0)
                {
                    mat_l=pow_gamma[norp][ser[1]][ser[2]];
                }
                else
                {
                    mat_l=pow_mgamma[norp][-ser[1]][ser[2]];
                }
            }
            else
            {
                if(ser[1]>0)
                {
                    mat_l=pow_q[norp][ser[1]-INT_MAX/2][ser[2]];
                }
                else
                {
                    mat_l=pow_qT[norp][-ser[1]-INT_MAX/2][ser[2]];
                }
                
            }
            double *mat_out=new double [sp_nljm_dim2[norp]];
            char side='r';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_r,&dim,mat_l,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(ser_r[0]>=mm_limit)
            {
                delete [] mat_r;
            }
            return mat_out;
        }
        int ser_r[5];
        ser_r[0]=4;
        memcpy(ser_r+1,ser+3,sizeof(int)*4);
        double *mat_r=mm_out(ser_r,norp);
        double *mat_l;
        if(abs(ser[1])<INT_MAX/2)
        {
            if(ser[1]>=0)
            {
                mat_l=pow_gamma[norp][ser[1]][ser[2]];
            }
            else
            {
                mat_l=pow_mgamma[norp][-ser[1]][ser[2]];
            }
        }
        else
        {
            if(ser[1]>0)
            {
                mat_l=pow_q[norp][ser[1]-INT_MAX/2][ser[2]];
            }
            else
            {
                mat_l=pow_qT[norp][-ser[1]-INT_MAX/2][ser[2]];
            }
            
        }
        double *mat_out=new double [sp_nljm_dim2[norp]];
        char opa='n';
        char opb='n';
        int dim=sp_nljm[norp].size();
        double alpha=1;
        double beta=0;
        // t3=currentcycles();
        dgemm_(&opa,&opb,&dim,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
        // t4=currentcycles();
        // t_sum+=t4-t3;
        if(ser_r[0]>=mm_limit)
        {
            delete [] mat_r;
        }
        return mat_out;
    }
    else if(ser[0]==8&&abs(ser[3])>=INT_MAX/2&&ser[3]==-ser[5]&&ser[4]==ser[6])
    {
        if(abs(ser[1])<INT_MAX/2&&ser[2]%2==0)
        {
            double *mat_l;
            if(ser[1]>=0)
            {
                mat_l=pow_gamma[norp][ser[1]][ser[2]];
            }
            else
            {
                mat_l=pow_mgamma[norp][-ser[1]][ser[2]];
            }
            int ser_r[7];
            ser_r[0]=6;
            memcpy(ser_r+1,ser+3,sizeof(int)*6);
            double *mat_r=mm_out(ser_r,norp);
            char side='l';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            double *mat_out=new double [sp_nljm_dim2[norp]];
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(ser_r[0]>=mm_limit)
            {
                delete [] mat_r;
            }
            return mat_out;
        }
        if(abs(ser[7])<INT_MAX/2&&ser[8]%2==0)
        {
            double *mat_r;
            if(ser[7]>=0)
            {
                mat_r=pow_gamma[norp][ser[7]][ser[8]];
            }
            else
            {
                mat_r=pow_mgamma[norp][-ser[7]][ser[8]];
            }
            int ser_l[7];
            ser_l[0]=6;
            memcpy(ser_l+1,ser+1,sizeof(int)*6);
            double *mat_l=mm_out(ser_l,norp);
            char side='r';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            double *mat_out=new double [sp_nljm_dim2[norp]];
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_r,&dim,mat_l,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(ser_l[0]>=mm_limit)
            {
                delete [] mat_l;
            }
            return mat_out;
        }
        double *mat_l;
        if(abs(ser[1])<INT_MAX/2)
        {
            if(ser[1]>=0)
            {
                mat_l=pow_gamma[norp][ser[1]][ser[2]];
            }
            else
            {
                mat_l=pow_mgamma[norp][-ser[1]][ser[2]];
            }
        }
        else
        {
            if(ser[1]>=0)
            {
                mat_l=pow_q[norp][ser[1]-INT_MAX/2][ser[2]];
            }
            else
            {
                mat_l=pow_qT[norp][-ser[1]-INT_MAX/2][ser[2]];
            }
        }
        int ser_r[7];
        ser_r[0]=6;
        memcpy(ser_r+1,ser+3,sizeof(int)*6);
        double *mat_r=mm_out(ser_r,norp);
        char opa='n';
        char opb='n';
        int dim=sp_nljm[norp].size();
        double alpha=1;
        double beta=0;
        double * mat_out=new double [sp_nljm_dim2[norp]];
        // t3=currentcycles();
        dgemm_(&opa,&opb,&dim,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
        // t4=currentcycles();
        // t_sum+=t4-t3;
        if(ser_r[0]>=mm_limit)
        {
            delete [] mat_r;
        }
        return mat_out;
    }
    else
    {
        int size_r=size/2/2*2;
        int size_l=size-size_r;
        if(abs(ser[size_l-1])>INT_MAX/2&&ser[size_l-1]==-ser[size_l+1]&&ser[size_l]==ser[size_l+2])
        {
            size_l-=2;
            size_r+=2;
        }
        int ser_l[size_l+1];
        int ser_r[size_r+1];
        ser_l[0]=size_l;
        ser_r[0]=size_r;
        memcpy(ser_l+1,ser+1,sizeof(int)*size_l);
        memcpy(ser_r+1,ser+size_l+1,sizeof(int)*size_r);
        double *mat_l=mm_out(ser_l,norp);
        double *mat_r=mm_out(ser_r,norp);
        if(abs(ser[1])>INT_MAX/2&&ser[1]==-ser[3]&&ser[2]==ser[4]&&size_l==4)
        {
            double *mat_out=new double [sp_nljm_dim2[norp]];
            char side='l';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(size_l>=mm_limit)
            {
                delete [] mat_l;
            }
            if(size_r>=mm_limit)
            {
                delete [] mat_r;
            }
            return mat_out;
        }
        if(abs(ser[size-3])>INT_MAX/2&&ser[size-3]==-ser[size-1]&&ser[size-2]==ser[size-1]&&size_r==4)
        {
            double *mat_out=new double [sp_nljm_dim2[norp]];
            char side='r';
            char uplo='u';
            int dim=sp_nljm[norp].size();
            double alpha=1;
            double beta=0;
            // t3=currentcycles();
            dsymm_(&side,&uplo,&dim,&dim,&alpha,mat_r,&dim,mat_l,&dim,&beta,mat_out,&dim);
            // t4=currentcycles();
            // t_sum+=t4-t3;
            if(size_l>=mm_limit)
            {
                delete [] mat_l;
            }
            if(size_r>=mm_limit)
            {
                delete [] mat_r;
            }
            return mat_out;
        }
        double *mat_out=new double [sp_nljm_dim2[norp]];
        char opa='n';
        char opb='n';
        int dim=sp_nljm[norp].size();
        double alpha=1;
        double beta=0;
        // t3=currentcycles();
        dgemm_(&opa,&opb,&dim,&dim,&dim,&alpha,mat_l,&dim,mat_r,&dim,&beta,mat_out,&dim);
        // t4=currentcycles();
        // t_sum+=t4-t3;
        if(size_l>=mm_limit)
        {
            delete [] mat_l;
        }
        if(size_r>=mm_limit)
        {
            delete [] mat_r;
        }
        return mat_out;
    }
}

struct node<int*,double *> *mm_root[2]={NULL,NULL};
int mm_cmp(int* &T1,int *&T2)
{
    if(T1[0]<T2[0])
    {
        return -1;
    }
    else if(T1[0]>T2[0])
    {
        return 1;
    }
    else
    {
        return memcmp(T1+1,T2+1,sizeof(int)*T1[0]);
    }
}

template<typename T_key,typename T_val> T_val mm_tree_search_insert(struct node<T_key,T_val> * &tree,T_key &k,int (*cmp)(T_key &,T_key&),T_val (*cal)(T_key &,int &),int &norp)
//k is elements, v is an element. k's cmp needs to be specified.
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
        tree->k=new int [k[0]+1];
        memcpy(tree->k,k,sizeof(int)*(k[0]+1));
		tree->l=NULL;
		tree->r=NULL;
        tree->v=cal(k,norp);//the calculation operation must be at the end of this tree search, or the key may be not evalued and in the cal operation, there may be some comparision to key.
		return tree->v;
	}
    // else
    // {
    //     node<T_key,T_val> *p=tree;
    //     node<T_key,T_val> *p_last=NULL;
    //     bool isL=true;
    //     while(p!=NULL)
    //     {
    //         if(cmp(k,tree->k)==0)
    //         {
    //             return p->v;
    //         }
    //         else if(cmp(k,tree->k)<0)
    //         {
    //             p_last=p;
    //             p=p->l;
    //             isL=true;
    //         }
    //         else
    //         {
    //             p_last=p;
    //             p=p->r;
    //             isL=false;
    //         }
    //     }
    //     if(isL)
    //     {
    //         p_last->l=(node<T_key,T_val>*) malloc(sizeof(struct node<T_key,T_val>));
    //         p_last->l->l=NULL;
    //         p_last->l->r=NULL;
    //         p_last->l->k=new int [k[0]+1];
    //         memcpy(p_last->l->k,k,sizeof(int)*(k[0]+1));
    //         p_last->l->v=cal(k,norp);//the calculation operation must be at the end of this tree search, or the key may be not evalued and in the cal operation, there may be some comparision to key.
    //         return p_last->l->v;
    //     }
    //     else
    //     {
    //         p_last->r=(node<T_key,T_val>*) malloc(sizeof(struct node<T_key,T_val>));
    //         p_last->r->l=NULL;
    //         p_last->r->r=NULL;
    //         p_last->r->k=new int [k[0]+1];
    //         memcpy(p_last->r->k,k,sizeof(int)*(k[0]+1));
    //         p_last->r->v=cal(k,norp);//the calculation operation must be at the end of this tree search, or the key may be not evalued and in the cal operation, there may be some comparision to key.
    //         return p_last->r->v;
    //     }
    // }
	if(cmp(k,tree->k)==0)
	{
		return tree->v;
	}
	else if(cmp(k,tree->k)<0)
	{
		return mm_tree_search_insert(tree->l,k,cmp,cal,norp);
	}
	else
	{
		return mm_tree_search_insert(tree->r,k,cmp,cal,norp);
	}
	
}
double *mm_out(int *ser,int &norp)
{
    if(ser[0]>=mm_limit)
    {
        return mm_cal(ser,norp);
    }
    else
    {
        return mm_tree_search_insert(mm_root[norp],ser,mm_cmp,mm_cal,norp);
    }
   
}