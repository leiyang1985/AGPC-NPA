void swap(int *a,int i,int j)
{
    a[i]^=a[j];
    a[j]^=a[i];
    a[i]^=a[j];
}

 

//将数组a中的下标i到下标j之间的所有元素逆序倒置
void reverse(int a[],int i,int j)
{
    for(; i<j; ++i,--j) 
    {
        swap(a,i,j);
    }
}
double O_cal(vector<unsigned char> &nt,int r[], int s[], int t,int nt_sum)
{
    if(nt[t]%2==1)
    {
        return -2;
    }
    else
    {
        return 2;
    }
}//这里我设所有的输入矩阵omega=（0，1，-1，0）二维的，那么omega^2=(1,0,0,1)为单位阵，由此得上述返回值
double piO(vector<unsigned char> &nt,int r[], int s[],int N)//该函数是对O（）连环乘积作连乘
{
    double res=1;
    int nt_sum=0;
    for(int t=0;t<nt.size();t++)
    {
        res*=O_cal(nt,r,s,t,nt_sum);
        nt_sum+=nt[t];
    }
    return res;
}
double s_com(vector<unsigned char> &nt,int r[], int s[],int N)//该函数对s(右矢)作全排列
{
    if(N<2)
    {
        //do something
        return piO(nt,r,s,N);
    }
    bool end=false;
    double res=0;
    while(true) 
    {
        res+=piO(nt,r,s,N);
        int i,j;
        //找到不符合趋势的元素的下标i
        for(i=N-2; i>=0; --i) 
        {
            if(s[i]<s[i+1]) 
            {
                break;
            }
            else if(i==0) 
            {
                return res;
            }
        }
        for(j=N-1; j>i; --j) 
        {
            if(s[j]>s[i])
            {
                break;
            }
        }
        swap(s,i,j);
        reverse(s,i+1,N-1);
    }
    // return res;
}


double r_com(vector<unsigned char> &nt,int r[],int N)//该函数对r(左矢)作全排列
{
    if(N<2)
    {
        int s[N];
        for(int i=0;i<N;i++)
        {
            s[i]=i;
        }
        return s_com(nt,r,s,N);
    }
    bool end=false;
    double res=0;
    while(true) 
    {
        int s[N];
        for(int i=0;i<N;i++)
        {
            s[i]=i;
        }
        res+=s_com(nt,r,s,N);//排列完全后，再调用s_com()进行s（右矢）全排列
        // for(int i=0;i<N;i++)
        // {
        //     cout<<r[i]<<' ';
        // }
        // cout<<endl;
        int i,j;
        //找到不符合趋势的元素的下标i
        for(i=N-2; i>=0; --i) 
        {
            if(r[i]<r[i+1]) 
            {
                break;
            }
            else if(i==0) 
            {
                return res;
            }
        }
        for(j=N-1; j>i; --j) 
        {
            if(r[j]>r[i])
            {
                break;
            }
        }
        swap(r,i,j);
        reverse(r,i+1,N-1);
    }
    // return res;
}


double com(int r[],int N)//该函数对r(左矢)作全排列
{
    if(N<2)
    {
        return 1;
    }
    bool end=false;
    double res=0;
    while(true) 
    {
        res++;
        //排列完全后，再调用s_com()进行s（右矢）全排列
        // for(int i=0;i<N;i++)
        // {
        //     cout<<r[i]<<' ';
        // }
        // cout<<endl;
        int i,j;
        //找到不符合趋势的元素的下标i
        for(i=N-2; i>=0; --i) 
        {
            if(r[i]<r[i+1]) 
            {
                break;
            }
            else if(i==0) 
            {
                return res;
            }
        }
        for(j=N-1; j>i; --j) 
        {
            if(r[j]>r[i])
            {
                break;
            }
        }
        swap(r,i,j);
        reverse(r,i+1,N-1);
    }
    //return res;
}

vector<vector<unsigned char>> nt_store;

void part(vector<unsigned char> &nt,int i_nt,int N_rest,int N)//对N_rest作划分，划分结果放入nt_store中
{
	int begin=1;
	if(i_nt>0)
	{
		begin=nt[i_nt-1];
	}
	for(int k=begin;k<=N_rest;k++)
	{
        if(i_nt>=nt.size())
        {
            nt.emplace_back(k);
        }
        else
        {
            nt[i_nt]=k;
        }
        
		N_rest-=k;
		if(N_rest>=k)
		{
			part(nt,i_nt+1,N_rest,N);
		}
		else if(N_rest==0)
		{
            vector<unsigned char> nt_temp(i_nt+1);
            copy(nt.begin(),nt.begin()+i_nt+1,nt_temp.begin());
            nt_store.emplace_back(nt_temp);
		}
        N_rest+=k;
	}
}



void ove_for_init(int N)//对N_rest=N作划分，划分结果放入nt_store中
{
    vector<unsigned char> nt;
    part(nt,0,N,N);
    int r[N];
    for(int i=0;i<N;i++)
    {
        r[i]=i;
    }
    double full_com_r=com(r,N);
    double m_sum=0;
    for(int i=0;i<nt_store.size();i++)
    {
        double m_delta=0;
        for(int k=0;k<nt_store[i].size();k++)
        {
            m_delta+=nt_store[i][k]*2-1;
        }
        int begin=0;
        int end=0;
        for(int l=0;l<nt_store[i].size();l++)
        {
            m_delta/=nt_store[i][l];//此处对应于(-1/2)^k 以及\Pi 1/n_t修正因子
            if(nt_store[i][l]!=nt_store[i][begin])//在序列中找与begin处不同的n_t
            {
                end=l-1;//找到后，设l-1为end，begin到end之间全为相同的n_t。
                m_delta/=fac_out(end-begin+1);//此处对应\Pi 1/R_nu!修正因子
                begin=l;
            }            
        }
        m_sum+=m_delta;
    }
    cout<<N<<' '<<m_sum<<' '<<full_com_r<<' '<<m_sum*full_com_r*full_com_r<<endl;
}

double ove_for_cal(int N)
{
    int r[N];
    double res=0;
    for(int t=0;t<nt_store.size();t++)//对每一个{n_t}划分进行历遍循环
    {
        for(int i=0;i<N;i++)
        {
            r[i]=i;
        }
        double delta=r_com(nt_store[t],r,N);//调用r_com()对r(左矢)进行全排列
        int begin=0;
        int end=0;
        for(int l=0;l<nt_store[t].size();l++)
        {
            delta*=-0.5/nt_store[t][l];//此处对应于(-1/2)^k 以及\Pi 1/n_t修正因子
            if(nt_store[t][l]!=nt_store[t][begin])//在序列中找与begin处不同的n_t
            {
                end=l-1;//找到后，设l-1为end，begin到end之间全为相同的n_t。
                delta/=fac_out(end-begin+1);//此处对应\Pi 1/R_nu!修正因子
                begin=l;
            }            
        }
        end=nt_store[t].size()-1;//此处明确最后一个n_t值的的R_nu归属
        delta/=fac_out(end-begin+1);//此处对应\Pi 1/R_nu!修正因子
        res+=delta;//此处对应于\sum_k求和
    }
    return res;
}
