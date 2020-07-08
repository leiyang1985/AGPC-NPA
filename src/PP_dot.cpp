int dot_cmp(int* &T1,int* &T2)
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
        return memcmp(T1+1,T2+1,sizeof(int)*(T1[0]));
    }
}
double dot_sum(node<int*,int> *&dot,int &norp)
{
	if(dot!=NULL)
	{
        double res=dot_sum(dot->l,norp);
        if(dot->v!=0)
        {
            res+=dot->v*tr_out(dot->k,norp);
        }
        delete [] dot->k;
		res+=dot_sum(dot->r,norp);
        free(dot);
        return res;
	}
    else
    {
        return 0;
    }
    
}

// template<typename T_key,typename T_val> double dot_sum(node<T_key,T_val> *&tree,int norp)//this is the loop verison, it has been proven unefficicy.
// {
// 	node<T_key,T_val> *p=tree;
// 	stack<node<T_key,T_val> *>* stack=NULL;
//     double res=0;
// 	while(p!=NULL||stack!=NULL)
// 	{
// 		while(p!=NULL)
// 		{
// 			push(stack,p);
// 			p=p->l;
// 		}
// 		if(stack!=NULL)
// 		{
// 			p=top(stack);
//             if(p->v!=0)
//             {
//                 res+=p->v*tr_out(p->k,norp);
//             }
//             delete [] p->k;
//             free(p);
// 			pop(stack);
// 			p=p->r;
// 		}
// 	}
//     return res;
// }

// void dot_del(node<int *,int> *&dot)
// {
//     if(dot!=NULL)
//     {
//         dot_del(dot->l);
//         dot_del(dot->r);
//         delete [] dot->k;
//         free(dot);
//     }
// }

// template<typename T_key,typename T_val> void dot_del(node<T_key,T_val> *&tree)
// {
// 	node<T_key,T_val> *p=tree;
// 	stack<node<T_key,T_val> *>* stack=NULL;
// 	while(p!=NULL||stack!=NULL)
// 	{
// 		while(p!=NULL)
// 		{
// 			push(stack,p);
// 			p=p->l;
// 		}
// 		if(stack!=NULL)
// 		{
// 			p=top(stack);
// 			delete [] p->k;
//             free(p);
// 			pop(stack);
// 			p=p->r;
// 		}
// 	}
// }

void dot_tree_search_add(struct node<int*,int> * &tree,int* k,int v,int (*cmp)(int* &,int* &))
{
    // if(tree==NULL)
    // {
    //     tree=(node<int*,int> *)malloc(sizeof(struct node<int*,int>));
	// 	tree->l=NULL;
	// 	tree->r=NULL;
	// 	tree->k=new int [k[0]+1];
    //     memcpy(tree->k,k,sizeof(int)*(k[0]+1));
    //     tree->v=v;
    // }
    // else
    // {
    //     node<int *,int> *p=tree;
    //     node<int *,int> *p_last=NULL;
    //     bool isL=true;
    //     while(p!=NULL)
    //     {
    //         if(cmp(k,p->k)==0)
    //         {
    //             p->v+=v;
    //             return;
    //         }
    //         else if(cmp(k,p->k)<0)
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
    //         p_last->l=(node<int*,int> *)malloc(sizeof(struct node<int *,int>));
    //         p_last->l->l=NULL;
    //         p_last->l->r=NULL;
    //         p_last->l->k=new int [k[0]+1];
    //         memcpy(p_last->l->k,k,sizeof(int)*(k[0]+1));
    //         p_last->l->v=v;
    //     }
    //     else
    //     {
    //         p_last->r=(node<int*,int> *)malloc(sizeof(struct node<int*,int>));
    //         p_last->r->l=NULL;
    //         p_last->r->r=NULL;
    //         p_last->r->k=new int [k[0]+1];
    //         memcpy(p_last->r->k,k,sizeof(int)*(k[0]+1));
    //         p_last->r->v=v;
    //     }        
    // }
    

	if(tree==NULL)
	{
		tree=(node<int*,int> *)malloc(sizeof(struct node<int*,int>));
		tree->l=NULL;
		tree->r=NULL;
		tree->k=new int [k[0]+1];
        memcpy(tree->k,k,sizeof(int)*(k[0]+1));
        tree->v=v;
	}
    else
    {
        if(cmp(k,tree->k)==0)
        {
            tree->v+=v;
        }
        else if(cmp(k,tree->k)<0)
        {
            dot_tree_search_add(tree->l,k,v,cmp);
        }
        else
        {
            dot_tree_search_add(tree->r,k,v,cmp);
        }
    }
}
double C_dot_sum(node<int*,int> *&dot,int &norp);
double PP_dot_cal(struct PQ_mat_struct & P1,struct PQ_mat_struct & P2,int norp)
{
    if(P1.M!=P2.M)
    {
        return 0;
    }
    if(P1.pi!=P2.pi)
    {
        return 0;
    }
    if(P1.J!=-1&&P2.J!=-1&&P1.J!=P2.J)
    {
        return 0;
    }
    bool isC=false;
    for(int i=1;i<=P1.ser[1][0];i+=2)
    {
        if(abs(P1.ser[1][i])>INT_MAX/4&&abs(P1.ser[1][i])<INT_MAX/2)
        {
            isC=true;
            break;
        }
    }
    for(int i=1;i<=P2.ser[1][0];i+=2)
    {
        if(abs(P2.ser[1][i])>INT_MAX/4&&abs(P2.ser[1][i])<INT_MAX/2)
        {
            isC=true;
            break;
        }
    }
    if(isC)
    {
        struct node<int*,int> * dot=NULL;
        for(int i=1;i<=P1.val[0];i++)
        {
            for(int j=1;j<=P2.val[0];j++)
            {
                int size1=P1.ser[i][0];
                int size2=P2.ser[j][0];
                int size=size1+size2;
                int val=P1.val[i]*P2.val[j];
                int ser[size+1];
                ser[0]=size;
                memcpy(ser+1,P1.ser[i]+1,sizeof(int)*(P1.ser[i][0]));
                if(P1.ser[i][P1.ser[i][0]-1]==P2.ser[j][1])
                {
                    ser[size1]+=P2.ser[j][2];
                    memcpy(ser+size1+1,P2.ser[j]+3,sizeof(int)*(size2-2));
                    ser[0]-=2;
                    size-=2;
                }
                else
                {
                    memcpy(ser+size1+1,P2.ser[j]+1,sizeof(int)*size2);
                }
                if(ser[1]==ser[size-1]&&size>2)
                {
                    ser[2]+=ser[size];
                    ser[0]-=2;
                    size-=2;
                }
                int C_index;
                for(int ii=1;ii<=size;ii+=2)
                {
                    if(abs(ser[ii])>INT_MAX/4&&abs(ser[ii])<INT_MAX/2)
                    {
                        C_index=ii;
                        break;
                    }
                }
                int ser_min[size+1];
                ser_min[0]=size;
                memcpy(ser_min+1,ser+C_index+2,sizeof(int)*(size-C_index-1));
                memcpy(ser_min+size-C_index,ser+1,sizeof(int)*(C_index+1));
                if(ser_min[1]>ser_min[size-1])
                {
                    int valT;
                    qT(val,ser_min,valT,ser);
                    dot_tree_search_add(dot,ser,valT,dot_cmp);
                }
                else
                {
                    dot_tree_search_add(dot,ser_min,val,dot_cmp);
                }
            }
        }
        return -0.5*C_dot_sum(dot,norp);
    }
    struct node<int*,int> * dot=NULL;
	for(int i=1;i<=P1.val[0];i++)
	{
		for(int j=1;j<=P2.val[0];j++)
		{
            int size1=P1.ser[i][0];
            int size2=P2.ser[j][0];
            int size=size1+size2;
			int val=P1.val[i]*P2.val[j];
            int ser[size+1];
            ser[0]=size;
            memcpy(ser+1,P1.ser[i]+1,sizeof(int)*(P1.ser[i][0]));
			if(P1.ser[i][P1.ser[i][0]-1]==P2.ser[j][1]&&(abs(P2.ser[j][1])<INT_MAX/4||abs(P2.ser[j][1])>INT_MAX/2))
			{
                ser[size1]+=P2.ser[j][2];
                memcpy(ser+size1+1,P2.ser[j]+3,sizeof(int)*(size2-2));
                ser[0]-=2;
                size-=2;
			}
			else
			{
                memcpy(ser+size1+1,P2.ser[j]+1,sizeof(int)*size2);
			}
			if(ser[1]==ser[size-1]&&size>2&&(abs(ser[1])<INT_MAX/4||abs(ser[1])>INT_MAX/2))
			{
                ser[2]+=ser[size];
                ser[0]-=2;
                size-=2;
			}
			if(size>2)
			{
                int ser_min[size+1];
                ser_min[0]=size;
                int k_begin=1;
                bool oddexist=false;
                for(int k=1;k<=size;k+=2)
                {
                    if(abs(ser[k])<INT_MAX/4&&ser[k+1]%2!=0)
                    {
                        //stop here**************
                        oddexist=true;
                        memcpy(ser_min+1,ser+k,sizeof(int)*(size-k+1));
                        memcpy(ser_min+size-k+2,ser+1,sizeof(int)*(k-1)); 
                        k_begin=k+2;
                        break;
                    }
                }
				int ser_con[size+1];
                ser_con[0]=size;
                if(oddexist)
                {
                    for(int k=k_begin;k<=size;k+=2)
                    {
                        if(abs(ser[k])<INT_MAX/2&&ser[k+1]%2!=0)
                        {
                            memcpy(ser_con+1,ser+k,sizeof(int)*(size-k+1));
                            memcpy(ser_con+size-k+2,ser+1,sizeof(int)*(k-1)); 
                            if(memcmp(ser_con+1,ser_min+1,sizeof(int)*size)<0)
                            {
                                memcpy(ser_min+1,ser_con+1,sizeof(int)*size);
                            }
                        }
                    }
                    int ser_T[size+1];
                    int val_T;
                    qT(val,ser,val_T,ser_T);
                    val=val_T;
                    memcpy(ser,ser_T,sizeof(int)*(size+1));
                    for(int k=1;k<=size;k+=2)
                    {
                        if(abs(ser[k])<INT_MAX/2&&ser[k+1]%2!=0)
                        {
                            memcpy(ser_con+1,ser+k,sizeof(int)*(size-k+1));
                            memcpy(ser_con+size-k+2,ser+1,sizeof(int)*(k-1));
                            if(memcmp(ser_con+1,ser_min+1,sizeof(int)*size)<0)
                            {
                                memcpy(ser_min+1,ser_con+1,sizeof(int)*size);
                            }
                        }
                    }
                }
                else
                {
                    for(int k=1;k<=size;k+=2)
                    {
                        if(abs(ser[k])<INT_MAX/2)
                        {
                            memcpy(ser_min+1,ser+k,sizeof(int)*(size-k+1));
                            memcpy(ser_min+size-k+2,ser+1,sizeof(int)*(k-1)); 
                            k_begin=k+2;
                            break;
                        }
                    }
                    for(int k=k_begin;k<=size;k+=2)
                    {
                        if(abs(ser[k])<INT_MAX/2)
                        {
                            memcpy(ser_con+1,ser+k,sizeof(int)*(size-k+1));
                            memcpy(ser_con+size-k+2,ser+1,sizeof(int)*(k-1));
                            if(memcmp(ser_con+1,ser_min+1,sizeof(int)*size)<0)
                            {
                                memcpy(ser_min+1,ser_con+1,sizeof(int)*size);
                            }
                        }
                    }
                    int ser_T[size+1];
                    int val_T;
                    qT(val,ser,val_T,ser_T);
                    val=val_T;
                    memcpy(ser,ser_T,sizeof(int)*(size+1));
                    for(int k=1;k<=size;k+=2)
                    {
                        if(abs(ser[k])<INT_MAX/2)
                        {
                            memcpy(ser_con+1,ser+k,sizeof(int)*(size-k+1));
                            memcpy(ser_con+size-k+2,ser+1,sizeof(int)*(k-1));
                            if(memcmp(ser_con+1,ser_min+1,sizeof(int)*size)<0)
                            {
                                memcpy(ser_min+1,ser_con+1,sizeof(int)*size);
                            }
                        }
                    }
                }
                memcpy(ser,ser_min,sizeof(int)*(size+1));
			}
            dot_tree_search_add(dot,ser,val,dot_cmp);
		}
	}
    return -0.5*dot_sum(dot,norp);
}