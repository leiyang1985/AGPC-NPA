int PQ_cmp(int * &T1,int *&T2)
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

void PQ_tree_op(node<int*,int> *&tree,int **&ser,int *&val,int &index_ser)
{
	// node<int *,int> *p=tree;
	// stack<node<int *,int> *>* stack=NULL;
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
	// 		if(p->v!=0)
	// 		{
	// 			ser[index_ser]=new int [p->k[0]+1];
	// 			memcpy(ser[index_ser],p->k,sizeof(int)*(p->k[0]+1));
	// 			val[index_ser]=p->v;
	// 			index_ser++;
	// 		}
	// 		delete [] p->k;
    //         free(p);
	// 		pop(stack);
	// 		p=p->r;
	// 	}
	// }
	
	
	if(tree!=NULL)
	{
		PQ_tree_op(tree->l,ser,val,index_ser);
		if(tree->v!=0)
		{
			ser[index_ser]=new int [tree->k[0]+1];
			memcpy(ser[index_ser],tree->k,sizeof(int)*(tree->k[0]+1));
			val[index_ser]=tree->v;
			index_ser++;
		}
		delete [] tree->k;
		PQ_tree_op(tree->r,ser,val,index_ser);
		free(tree);
	}
}
template<typename T_key,typename T_val> void PQ_tree_search_add(struct node<T_key,T_val> * &tree,T_key k,T_val &v,int (*cmp)(T_key &,T_key&))
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
        tree->k=new int [k[0]+1];
		memcpy(tree->k,k,sizeof(int)*(k[0]+1));
        tree->v=v;
	}
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
	else
    {
        if(cmp(k,tree->k)==0)
        {
            tree->v+=v;
        }
        else if(cmp(k,tree->k)<0)
        {
            PQ_tree_search_add(tree->l,k,v,cmp);
        }
        else
        {
            PQ_tree_search_add(tree->r,k,v,cmp);
        }
    }
}
PQ_mat_struct PQ_con_cal(PQ_mat_struct &P,PQ_mat_struct &Q,int &coe)
{
	PQ_mat_struct temp;
	struct node<int *,int> *PQ_root=NULL;
	
	temp.M=P.M-Q.M;
	temp.pi=(P.pi+Q.pi)%2;
	if(P.J==0||Q.J==0)
	{
		temp.J=P.J+Q.J;
	}
	else
	{
		temp.J=-1;
	}
	temp.PQ_Jpi_index=-1;
	for(int i=1;i<=P.val[0];i++)
	{
		for(int j=1;j<=Q.val[0];j++)
		{
			int size_P=P.ser[i][0];
			int size_Q=Q.ser[j][0];
			int size=P.ser[i][0]+Q.ser[j][0];
			int ser[size+1];
			ser[0]=size;
			int val=P.val[i]*Q.val[j];
			memcpy(ser+1,P.ser[i]+1,sizeof(int)*size_P);
			if(P.ser[i][size_P-1]==Q.ser[j][1])
			{
				ser[0]-=2;
				size-=2;
				ser[size_P]+=Q.ser[j][2];
				memcpy(ser+size_P+1,Q.ser[j]+3,sizeof(int)*(size_Q-2));
			}
			else
			{
				memcpy(ser+size_P+1,Q.ser[j]+1,sizeof(int)*size_Q);
			}
			PQ_tree_search_add(PQ_root,ser,val,PQ_cmp);
			int ser_T[size+1];
			int val_T;
			qT(val,ser,val_T,ser_T);
			val_T*=-1;
			PQ_tree_search_add(PQ_root,ser_T,val_T,PQ_cmp);
		}
	}
	int num_ser=2*P.val[0]*Q.val[0];
	temp.val=new int [num_ser+1];
	temp.ser=new int *[num_ser+1];
	temp.ser[0]=NULL;
	int index_ser=1;
	PQ_tree_op(PQ_root,temp.ser,temp.val,index_ser);
	temp.val[0]=index_ser-1;
	if(temp.val[0]!=0)
	{
		coe=abs(temp.val[1]);
		for(int i=2;i<=temp.val[0];i++)
		{
			coe=gcb(coe,temp.val[i]);
		}
		if(temp.val[1]<0)
		{
			coe=-coe;
		}
		for(int i=1;i<=temp.val[0];i++)
		{
			temp.val[i]/=coe;
		}
	}
	else
	{
		coe=0;
	}
	return temp;
}