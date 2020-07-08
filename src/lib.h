template<typename T> void swap(T *a,T *b)
{
	T temp=*a;
	*a=*b;
	*b=temp;
}
inline int pack_k(int i,int j)
{
	return i+(j+1)*j/2;
}
typedef unsigned long long cycles_t; 
inline cycles_t currentcycles() 
{ 
     __u32 lo,hi;

        __asm__ __volatile__
        (
         "rdtsc":"=a"(lo),"=d"(hi)
        );
        return (__u64)hi<<32|lo;
} 
int cmp(const void *p1,const void *p2)
{
	int i1=*((int *)p1);
	int i2=*((int *)p2);
	if(i1==i2)
	{
		return 0;
	}
	else if(i1<i2)
	{
		return -1;
	}
	else
	{
		return 1;
	}
}

extern "C" double ddot_(int *N,double *DX,int *INCX,double *DY,int *INCY);
extern "C" double dnrm2_(int *N,double *x,int *ins);
extern "C" void dcopy_(int *N,double *from,int *ins_from,double *to,int *ins_to);
extern "C" void daxpy_(int *N,double *alpha,double *x,int *incx,double *y,int *incy);
extern "C" void dscal_(int *N,double *alpha,double *v_inout,int *ins);
//extern "C" void dsyev_(char *a,char *b,int *n,double *A,int *nn,double *e,double *work,int *lwork,int *info);
extern "C" void dsymm_(char *side,char *uorl,int *M,int *N,double *alpha,double *A,int *lda,double *B,int *ldb,double *beta,double *C,int *ldc);
extern "C" void dsymv_(char *uplo,int *N, double *alpha,double *A,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
extern "C" void dgemm_(char *TRANSA,char *TRANSB,int *M,int *N,int *K,double *ALPHA,double *A,int *LDA,double *B,int *LDB,double *BETA,double *C,int *LDC);
extern "C" void dgemv_(char *TRANS, int *M, int *N, double *ALPHA, double *A, int *LDA, double *X, int *INCX,double *BETA, double *Y,int *INCY);
//extern "C" void dgesvd_(char * JOBU,char * JOBVT,int * M,int * N, double *A, int *LDA, double *S, double *U, int *LDU, double *VT,int *  LDVT, double *WORK, int *LWORK, int *INFO );
extern "C" double * TRANSPOSE_(double *array);
//extern "C" void mkl_dimatcopy(const char ordering, const char trans, size_t rows, size_t cols, const double alpha, double * AB, size_t lda, size_t ldb);
template<typename T_key,typename T_val> struct node
{
	T_key k;
	T_val v;
	struct node<T_key,T_val> *l;
	struct node<T_key,T_val> *r;
};

template<typename T_key,typename T_val> T_val tree_search_insert(struct node<T_key,T_val> * &tree,T_key &k,int (*cmp)(T_key &,T_key&,void *),void *par_cmp,void (*evalue)(T_key&,T_key&,void *),void *par_evalue,void (*cal)(T_val &,T_key &,void *),void *par_cal,void(*del)(T_key &,void *),void *par_del)
//most general case
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
        memset(tree,0,sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
		cal(tree->v,k,par_cal);
		evalue(tree->k,k,par_evalue);
		return tree->v;
	}
	if(cmp(k,tree->k,par_cmp)==0)
	{
		del(k,par_del);
		return tree->v;
	}
	if(cmp(k,tree->k,par_cmp)<0)
	{
		return tree_search_insert(tree->l,k,cmp,par_cmp,evalue,par_evalue,cal,par_cal,del,par_del);
	}
	if(cmp(k,tree->k,par_cmp)>0)
	{
		return tree_search_insert(tree->r,k,cmp,par_cmp,evalue,par_evalue,cal,par_cal,del,par_del);
	}
	
}
template<typename T_key,typename T_val> T_val tree_search_insert(struct node<T_key,T_val> * &tree,T_key &k,int (*cmp)(T_key &,T_key&),T_val (*cal)(T_key &))
//k is a seiral of elements, v is an element, and k need to be deleted if found. k's cmp needs to be specified.
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
        memset(tree,0,sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
        tree->v=cal(k);
		tree->k=k;
		return tree->v;
	}
	if(cmp(k,tree->k)==0)
	{
		delete [] k;
		return tree->v;
	}
	if(cmp(k,tree->k)<0)
	{
		return tree_search_insert(tree->l,k,cmp,cal);
	}
	if(cmp(k,tree->k)>0)
	{
		return tree_search_insert(tree->r,k,cmp,cal);
	}
	
}
template<typename T_key,typename T_val> T_val tree_search_insert(struct node<T_key,T_val> * &tree,T_key &k,T_val (*cal)(T_key&))
//this is for k is an element, and comparabe by <>= and v is an element
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
        memset(tree,0,sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
        tree->v=cal(k);
        tree->k=k;
		return tree->v;
	}
	if(k==tree->k)
	{
		return tree->v;
	}
	if(k<tree->k)
	{
		return tree_search_insert(tree->l,k,cal);
	}
	if(k>tree->k)
	{
		return tree_search_insert(tree->r,k,cal);
	}
	
}
template<typename T_key,typename T_val> T_val tree_search_insert(struct node<T_key,T_val> * &tree,T_key &k,int (*cmp)(T_key &,T_key&),T_val (*cal)(T_key &,void *),void *par_cal)
//k is elements, v is an element. k's cmp needs to be specified.
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
        memset(tree,0,sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
        tree->v=cal(k,par_cal);
		tree->k=k;
		return tree->v;
	}
	if(cmp(k,tree->k)==0)
	{
		return tree->v;
	}
	if(cmp(k,tree->k)<0)
	{
		return tree_search_insert(tree->l,k,cmp,cal,par_cal);
	}
	if(cmp(k,tree->k)>0)
	{
		return tree_search_insert(tree->r,k,cmp,cal,par_cal);
	}
	
}
template<typename T_key,typename T_val> void tree_search_add(struct node<T_key,T_val> * &tree,T_key &k,T_val &v,int (*cmp)(T_key &,T_key&,void *),void *par_cmp,void (*evalue_k)(T_key&,T_key&,T_val &,T_val &,void *),void *par_evalue_k,void (*evalue_v)(T_key&,T_key&,T_val &,T_val &,void *),void *par_evalue_v,void (*add)(T_val &,T_val &,void *),void *par_add,void(*del)(T_key &,void *),void *par_del)
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
        memset(tree,0,sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
		evalue_k(tree->k,k,par_evalue_k);
		evalue_v(tree->v,v,par_evalue_v);
	}
    else
    {
        if(cmp(k,tree->k,par_cmp)==0)
        {
            add(tree->v,v,par_add);
        }
        if(cmp(k,tree->k,par_cmp)<0)
        {
            tree_search_add(tree->l,k,v,cmp,par_cmp,evalue_k,par_evalue_k,evalue_v,par_evalue_v,add,par_add,del,par_del);
        }
        if(cmp(k,tree->k,par_cmp)>0)
        {
            tree_search_add(tree->r,k,v,cmp,par_cmp,evalue_k,par_evalue_k,evalue_v,par_evalue_v,add,par_add,del,par_del);
        }
    }
}
template<typename T>int ser_cmp(vector<T> &T1,vector<T> &T2)
{
    if(T1.size()<T2.size())
    {
        return -1;
    }
    else if(T1.size()>T2.size())
    {
        return 1;
    }
    else
    {
        return(&(T1[0]),&(T2[0]),sizeof(T)*T1.size());
    }
}
template<typename T_key,typename T_val> void tree_search_add(struct node<T_key,T_val> * &tree,T_key &k,T_val v,int (*cmp)(T_key &,T_key&))
{
	if(tree==NULL)
	{
		tree=(node<T_key,T_val> *)malloc(sizeof(struct node<T_key,T_val>));
        memset(tree,0,sizeof(struct node<T_key,T_val>));
		tree->l=NULL;
		tree->r=NULL;
        tree->k=k;
        tree->v=v;
	}
    else
    {
        if(cmp(k,tree->k)==0)
        {
            tree->v+=v;
        }
        if(cmp(k,tree->k)<0)
        {
            tree_search_add(tree->l,k,v,cmp);
        }
        if(cmp(k,tree->k)>0)
        {
            tree_search_add(tree->r,k,v,cmp);
        }
    }
}

template<typename T_key,typename T_val> void tree_op(node<T_key,T_val> *&tree,void(*op)(node<T_key,T_val> *&,void *),void *par_op)
{
	if(tree!=NULL)
	{
		tree_op(tree->l,op,par_op);
		op(tree,par_op);
		tree_op(tree->r,op,par_op);
	}
}

template<typename T_key,typename T_val> inline void key_val_init(T_key* &key,T_val* &val,size_t &dim,size_t &dim_max)
{
	dim_max=256;
	dim=0;
	key=new T_key [dim_max];
	val=new T_val [dim_max];
}
template <typename T> bool array_search(T* array,T target,size_t & index,size_t dim, int (*cmp)(T&,T&,void *),void *cmp_par)
{
	if(dim<=0)
	{
		index=0;
		return(false);
	}
	size_t begin,end,mid;
	begin=0;
	end=dim-1;
	if(cmp(array[begin],target,NULL)>0)
	{
		index=begin;
		return false;
	}
	else if(cmp(array[begin],target,NULL)==0)
	{
		index=begin;
		return true;
	}
	if(cmp(array[end],target,NULL)<0)
	{
		index=end+1;
		return false;
	}
	else if(cmp(array[end],target,NULL)==0)
	{
		index=end;
		return true;
	}
	while(begin+2<=end)
	{
		mid=(end+begin)/2;
		if(cmp(array[mid],target,NULL)<0)
		{
			begin=mid;
			continue;
		}
		else if(cmp(array[mid],target,NULL)>0)
		{
			end=mid;
			continue;
		}
		else
		{
			index=mid;
			return true;
		}		
	}
	index=end;
	return false;
}
template <typename T> bool array_search(T*  array,T target,size_t & index,size_t dim, int (*cmp)(T&,T&))
//cmp does not need par.
{
	if(dim<=0)
	{
		index=0;
		return(false);
	}
	size_t begin,end,mid;
	begin=0;
	end=dim-1;
	if(cmp(array[begin],target)>0)
	{
		index=begin;
		return false;
	}
	else if(cmp(array[begin],target)==0)
	{
		index=begin;
		return true;
	}
	if(cmp(array[end],target)<0)
	{
		index=end+1;
		return false;
	}
	else if(cmp(array[end],target)==0)
	{
		index=end;
		return true;
	}
	while(begin+2<=end)
	{
		mid=(end+begin)/2;
		if(cmp(array[mid],target)<0)
		{
			begin=mid;
			continue;
		}
		else if(cmp(array[mid],target)>0)
		{
			end=mid;
			continue;
		}
		else
		{
			index=mid;
			return true;
		}		
	}
	index=end;
	return false;
}

template<typename T> int cmp(vector<T >&T1, vector<T >&T2)
{
	if(T1.size()<T2.size())
	{
		return -1;
	}
	else if(T1.size()>T2.size())
	{
		return 1;
	}
	else
	{
		return memcmp(&(T1[0]),&(T2[0]),sizeof(T )*T1.size());
	}
	
}

template <typename T> bool array_search(vector<T>  array,T target,size_t & index,size_t dim, int (*cmp)(T&,T&))
//cmp does not need par.
{
	if(dim<=0)
	{
		index=0;
		return(false);
	}
	size_t begin,end,mid;
	begin=0;
	end=dim-1;
	if(cmp(array[begin],target)>0)
	{
		index=begin;
		return false;
	}
	else if(cmp(array[begin],target)==0)
	{
		index=begin;
		return true;
	}
	if(cmp(array[end],target)<0)
	{
		index=end+1;
		return false;
	}
	else if(cmp(array[end],target)==0)
	{
		index=end;
		return true;
	}
	while(begin+2<=end)
	{
		mid=(end+begin)/2;
		if(cmp(array[mid],target)<0)
		{
			begin=mid;
			continue;
		}
		else if(cmp(array[mid],target)>0)
		{
			end=mid;
			continue;
		}
		else
		{
			index=mid;
			return true;
		}		
	}
	index=end;
	return false;
}

template <typename T_key,typename T_val> inline void key_val_double(T_key *&key,T_val *&val,size_t &dim)
{
	T_key key_temp[dim];
	T_val val_temp[dim];
	memcpy(key_temp,key,sizeof(T_key)*dim);
	memcpy(val_temp,val,sizeof(T_val)*dim);
	delete [] key;
	delete [] val;
	key=new T_key [dim*2];
	val=new T_val [dim*2];
	memcpy(key,key_temp,sizeof(T_key)*dim);	
	memcpy(val,val_temp,sizeof(T_val)*dim);
	dim*=2;
}

template <typename T_key,typename T_val> inline T_val key_val_cal_insert(T_key *&key,T_val *&val,T_key &key_insert,size_t &index,size_t &dim,size_t &dim_max,void (*evalue)(T_key&,T_key&,void *),void *par_evalue,void (*cal)(T_val &,T_key &,void *),void *par_cal)
{
	if(dim>index)
	{
		memmove(key+index+1,key+index,sizeof(T_key)*(dim-index));
		memmove(val+index+1,val+index,sizeof(T_val)*(dim-index));
	}
	evalue(key[index],key_insert,par_evalue);
	cal(val[index],key_insert,par_cal);
	dim++;
	if(dim>=dim_max)
	{
		key_val_double(key,val,dim_max);
	}
	return val[index];
}
template<typename T_key,typename T_val> inline T_val val_out(T_key *&key,T_val *&val,T_key &key_input,size_t &dim,size_t &dim_max,int (*cmp)(T_key &, T_key&,void *),void *cmp_par,void (*evalue)(T_key&,T_key&,void *),void *par_evalue,void (*cal)(T_val &,T_key &,void *),void *par_cal,void (*del)(T_key &,void *),void *del_par)
{
	size_t index;
	if(array_search(key,key_input,index,dim,cmp,cmp_par))
	{
		del(key_input,del_par);
		return val[index];
	}
	else
	{
		return key_val_cal_insert(key,val,key_input,index,dim,dim_max,evalue,par_evalue,cal,par_cal);
	}
}
unsigned int gcd(unsigned int a,unsigned int b)
{
    while(b^=a^=b^=a%=b);
    return a;
}

unsigned int gcb(unsigned int numberA, unsigned int numberB)
//this function is to obtain greatest common divisor of two unsingned ints. So make these two nozero first.
{
    if(numberA == numberB)
        return numberA;
    if(numberA < numberB)
        return gcd(numberB, numberA);
    else
    {
        //和1做按位与运算，判断奇偶
        if(!numberA&1 && !numberB&1)  //都是偶数
            return gcd(numberA>>1, numberB>>1) << 1;
        else if(!numberA&1 && numberB&1)
            return gcd(numberA>>1, numberB);
        else if(numberA&1 && !numberB&1)
            return gcd(numberA, numberB>>1);
        else
            return gcd(numberB, numberA - numberB);
    }
}
void lapack(double **h,double **vec,double *e,int dim)
{
	if(dim==1)
	{
		vec[0][0]=1;
		e[0]=h[0][0];
	}
	else
	{
		char a='v';
		char b='L';
		int N1;
		int NN;
		int lwork;
		int info;
		lwork=dim*3;
		double *work=new double [dim*3];
		N1=dim;
		NN=N1;
		info=0;
		double *h_temp=new double [dim*dim];
		int i,j;
		for(i=0;i<=dim-1;i++)
		{
			for(j=0;j<=dim-1;j++)
			{
				h_temp[i*dim+j]=h[i][j];
			}
		}
		dsyev_(&a,&b,&N1,h_temp,&NN,e,work,&lwork,&info);
		for(i=0;i<=dim-1;i++)
		{
			for(j=0;j<=dim-1;j++)
			{
				vec[i][j]=h_temp[i*dim+j];//each column in the h_temp in fortran style is a vec
			}
		}
		delete [] h_temp;
		delete [] work;
	}
}

template<typename T> struct stack
{
	T data;
	struct stack *next;
};

template<typename T> void push(stack<T> *&stack_current,T data)
{
	stack<T> *p=(stack<T> *)malloc(sizeof(stack<T>));
	p->data=data;//formally
	p->next=stack_current;
	stack_current=p;
}
template<typename T> void pop(stack<T> *&stack_current)
{
	//free(stack_crrent->data);//formally
	//delete [] stack_current->data;
	stack<T> *p_last=stack_current;
	stack_current=stack_current->next;
	free(p_last);
}
template<typename T> T top(stack<T> *&stack_current)
{
	return stack_current->data;
}
void geev_lapack(double *A,double *B,int dim,double *e)
{
	int itype=1;
	char uplo='U';
	char jobz='V';
	double *work=new double [6*dim];
	int lwork=6*dim;
	int info;
	dsygv_(&itype,&jobz,&uplo,&dim,A,&dim,B,&dim,e,work,&lwork,&info);//带ovelape 的对角化
	delete [] work;
	// work=new double [lwork];
	// dsygv_(&itype,&jobz,&uplo,&dim,A,&dim,B,&dim,e,work,&lwork,&info);
	// delete [] work;
	cout<<info<<endl;
}
template<typename T_key,typename T_val> void tree_op_loop(node<T_key,T_val> *&tree,void(*op)(node<T_key,T_val> *&,void *),void *par_op)
{
	node<T_key,T_val> *p=tree;
	stack<node<T_key,T_val> *>* stack=NULL;
	while(p!=NULL||stack!=NULL)
	{
		while(p!=NULL)
		{
			push(stack,p);
			p=p->l;
		}
		if(stack!=NULL)
		{
			p=top(stack);
			op(p,par_op);
			pop(stack);
			p=p->r;
		}
	}
}

/*
int svd_demo_lapack()
{
    int m=5;
    int n=10;
    int min,max;
    if(m<n)
    {
        min=m;
        max=n;
    }
    else
    {
        min=n;
        max=m;
    }
    if(min<2)
    {
        cout<<"you matrix has no svd possibility"<<endl;
        return 0;
    }
    double *a=new double [m*n];
    for(int i=0;i<m*n;i++)
    {
        a[i]=gauss(0,1);
    }
    double *a_bak=new double [m*n];
    memcpy(a_bak,a,sizeof(double)*m*n);
    int lda=m;
    int lwork=-1;
    char jobu='S';
    char jobvt='S';
    double *s=new double [min];
    double *u=new double [m*min];
    int ldu=m;
    double *vt=new double [min*n];
    int ldvt=min;
    double *work=new double [1];
    lwork=-1;
    int info;
    dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);
    lwork=work[0];
    delete [] work;
    work=new double [lwork];
    dgesvd_(&jobu,&jobvt,&m,&n,a,&lda,s,u,&ldu,vt,&ldvt,work,&lwork,&info);
    double a_final[m][n];
    cout<<"s=";
    for(int imin=0;imin<min;imin++)
    {
        cout<<s[imin]<<' ';
    }
    cout<<endl;
    cout<<"u="<<endl;
    for(int im=0;im<m;im++)
    {
        for(int imin=0;imin<min;imin++)
        {
            cout<<u[im+imin*m]<<' ';
        }
        cout<<endl;
    }
    cout<<"vt="<<endl;
    for(int imin=0;imin<min;imin++)
    {
        for(int in=0;in<n;in++)
        {
            cout<<vt[imin+in*min]<<' ';
        }
        cout<<endl;
    }
    mkl_dimatcopy('c','t',min,n,1,vt,min,n);//we need tranvers vt to make use to easily unserstand
    cout<<"v="<<endl;
    for(int in=0;in<n;in++)
    {
        for(int imin=0;imin<min;imin++)
        {
            cout<<vt[in+imin*n]<<' ';
        }
        cout<<endl;
    }
    cout<<"a_final="<<endl;
    for(int im=0;im<m;im++)
    {
        for(int in=0;in<n;in++)
        {
            a_final[im][in]=0;
            for(int k=0;k<min;k++)
            {
                a_final[im][in]+=s[k]*u[k*m+im]*vt[k*n+in];//each column of u or vt is a vec.
            }
            cout<<a_final[im][in]<<' ';
        }
        cout<<endl;
    }
    cout<<"a="<<endl;
     for(int im=0;im<m;im++)
    {
        for(int in=0;in<n;in++)
        {
            cout<<a_bak[im+in*m]<<' ';
        }
        cout<<endl;
    }
    delete [] work;
    delete [] a;
    delete [] a_bak;
    delete [] s;
    delete [] vt;
    delete [] u;
    return 0;
}
*/