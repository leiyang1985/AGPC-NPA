#define tr_limit 12

struct node<int *,double> *tr_root[2]={NULL,NULL};
double tr_cal(int *&ser,int &norp)
{	
    int size=ser[0];
	if(size==4)
	{
        double *mat_r;
        double *mat_l;
        if(abs(ser[1])<INT_MAX/4)
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
        if(abs(ser[3])<INT_MAX/4)
        {
            if(ser[3]>0)
            {
                mat_r=pow_gamma[norp][ser[3]][ser[4]];
            }
            else
            {
                mat_r=pow_mgamma[norp][-ser[3]][ser[4]];
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
        int dim=sp_nljm_dim2[norp];
        int inc=1;
        // t3=currentcycles();
        double res=ddot_(&dim,mat_l,&inc,mat_r,&inc);
        // t4=currentcycles();
        // t_sum+=t4-t3;
        if(ser[2]%2!=0)
        {
            res=-res;
        }
        return res;
    }
    else
    {
        double *mat_r;
        double *mat_l;
        if(ser[1]>=0)
        {
            mat_l=pow_gamma[norp][ser[1]][ser[2]];
        }
        else
        {
            mat_l=pow_mgamma[norp][-ser[1]][ser[2]];
        }
        int ser_low[size-1];
        ser_low[0]=size-2;
        memcpy(ser_low+1,ser+3,sizeof(int)*(size-2));        
        mat_r=mm_out(ser_low,norp);
        int dim=sp_nljm_dim2[norp];
        int inc=1;
        // t3=currentcycles();
        double res=ddot_(&dim,mat_l,&inc,mat_r,&inc);
        // t4=currentcycles();
        // t_sum+=t4-t3;
        if(ser_low[0]>=mm_limit)
        {
            delete [] mat_r;
        }
        if(ser[2]%2!=0)
        {
            res=-res;
        }
        return res;
	}
}
double C_tr_cal(int *&ser,int &norp)
{	
    int size=ser[0];
    if(size==4)
    {
        if(ser[1]==ser[3])
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else if(size==6)
    {
        double *mat;
        if(abs(ser[3])<INT_MAX/4)
        {
            if(ser[3]>0)
            {
                mat=pow_gamma[norp][ser[3]][ser[4]];
            }
            else
            {
                mat=pow_mgamma[norp][-ser[3]][ser[4]];
            }
        }
        else
        {
            
            if(ser[3]>0)
            {
                mat=pow_q[norp][ser[3]-INT_MAX/2][ser[4]];  
            }
            else
            {
                mat=pow_qT[norp][-ser[3]-INT_MAX/2][ser[4]];  
            }    
        }
        int ii,jj;
        if(ser[1]>0)
        {
            ii=ser[1]-INT_MAX/4-1;
        }
        else
        {
            ii=sp_nljm[norp].size()+ser[1]+INT_MAX/4;
        }
        if(ser[5]>0)
        {
            jj=ser[5]-INT_MAX/4-1;
        }
        else
        {
            jj=sp_nljm[norp].size()+ser[5]+INT_MAX/4;
        }
        return mat[ii+jj*sp_nljm[norp].size()];
    }
    else
    {
        double *mat;
        int ser_low[size-3];
        ser_low[0]=size-4;
        memcpy(ser_low+1,ser+3,sizeof(int)*(size-4));        
        mat=mm_out(ser_low,norp);
        int ii,jj;
        if(ser[1]>0)
        {
            ii=ser[1]-INT_MAX/4-1;
        }
        else
        {
            ii=sp_nljm[norp].size()+ser[1]+INT_MAX/4;
        }
        if(ser[size-1]>0)
        {
            jj=ser[size-1]-INT_MAX/4-1;
        }
        else
        {
            jj=sp_nljm[norp].size()+ser[size-1]+INT_MAX/4;
        }
        double res=mat[ii+jj*sp_nljm[norp].size()];
        if(ser_low[0]>=mm_limit)
        {
            delete [] mat;
        }
        return res;
	}
}
int tr_cmp(int* &T1,int * &T2)
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
template<typename T_key,typename T_val> T_val tr_tree_search_insert(struct node<T_key,T_val> * &tree,T_key &k,int (*cmp)(T_key &,T_key&),T_val (*cal)(T_key &,int&),int &norp)
//k is elements, v is an element. k's cmp needs to be specified.
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
        tree->k=new int [k[0]+1];
        memcpy(tree->k,k,sizeof(int)*(k[0]+1));
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
		return tr_tree_search_insert(tree->l,k,cmp,cal,norp);
	}
	else
	{
		return tr_tree_search_insert(tree->r,k,cmp,cal,norp);
	}
	
}
double tr_out(int* &ser,int norp)
{
    int size=ser[0];
    if(size==2)
	{
		if(abs(ser[1])<INT_MAX/4)
		{
			if(ser[1]>=0)
			{
                return tr_pow_gamma[norp][ser[1]][ser[2]];
                //t3=currentcycles();
				// double res=0;
				// for(int i=0;i<sp_nljm[norp].size();i++)
				// {
				// 	res+=pow_gamma[norp][ser[0]][ser[1]][i+i*sp_nljm[norp].size()];
				// }
                //t4=currentcycles();
               //t_sum+=t4-t3;
				// return res;
			}
			else
			{
                return tr_pow_mgamma[norp][-ser[1]][ser[2]];
                //t3=currentcycles();
				// double res=0;
				// for(int i=0;i<sp_nljm[norp].size();i++)
				// {
				// 	res+=pow_mgamma[norp][-ser[0]][ser[1]][i+i*sp_nljm[norp].size()];
				// }
                //t4=currentcycles();
               //t_sum+=t4-t3;
				// return res;
			}			
		}
		else
		{
            //t3=currentcycles();
			//int index=abs(ser[0])-INT_MAX/2;
            return tr_pow_q[norp][abs(ser[1])-INT_MAX/2][ser[2]];
			// double res=0;
			// for(int i=0;i<sp_nljm[norp].size();i++)
			// {
			// 	res+=pow_q[norp][index][ser[1]][i+i*sp_nljm[norp].size()];
			// }
            //t4=currentcycles();
           //t_sum+=t4-t3;
			// return res;
		}
	}
    return tr_cal(ser,norp);
    // if(size>=tr_limit)
    // {
    //     return tr_cal(ser,norp);
    // }
    // else
    // {
    //     return tr_tree_search_insert(tr_root[norp],ser,tr_cmp,tr_cal,norp);
    // }
}