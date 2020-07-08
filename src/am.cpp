#define PI 3.14159265358979323846

struct node<int,double> *fac_root=NULL;

double fac_out(int m);
inline double fac_cal(int &m)
{
	if(m<=1)
	{
		return 1;
	}
	else
	{
		return m*fac_out(m-1);
	}
}

inline double fac_out(int m)
{
	return tree_search_insert(fac_root,m,fac_cal);
}


struct node<int,double> *fac2_root=NULL;

double fac2_out(int m);
inline double fac2_cal(int &m)
{
	if(m<=1)
	{
		return 1;
	}
	else
	{
		return m*fac2_out(m-2);
	}
}

inline double fac2_out(int m)
{
	return tree_search_insert(fac2_root,m,fac2_cal);
}

inline int sgn(int i_phase)
{
	if(i_phase%2==0)
	{
		return 1;
	}
	else
	{
		return -1;
	}
}


inline bool threej_zero(int *j)
//using the order of {j1,j2,j3,m1,m2,m3}
{
	if(j[0]<0||j[1]<0||j[2]<0)
	{
		return true;
	}
	if(j[3]<-j[0]||j[3]>j[0])
	{
		return true;
	}
	if(j[4]<-j[1]||j[4]>j[1])
	{
		return true;
	}
	if(j[5]<-j[2]||j[5]>j[5])
	{
		return true;
	}
	if(j[3]+j[4]+j[5]!=0)
	{
		return true;
	}
	if(j[2]<abs(j[0]-j[1])||j[2]>j[0]+j[1])
	{
		return true;
	}
	if((j[0]+j[1]+j[2])%2!=0)
	{
		return true;
	}
	return false;
}


double cg_cal(int J1, int m1, int J2,int m2, int J3, int m3)
{
	int v_min,v_max;
	v_max=(J1+J2-J3);
	if((J1-m1)<v_max) v_max=J1-m1;
	if((J2+m2)<v_max) v_max=J2+m2;
	v_min=0;
	if((J1+m2-J3)>v_min) v_min=J1+m2-J3;
	if((J2-J3-m1)>v_min) v_min=J2-J3-m1;
	double sum=0;
	double delta;
	for(int v=v_min; v<=v_max; v+=2)
	{
		
		delta=fac_out(v/2)*fac_out((J1+J2-J3-v)/2)*fac_out((J1-m1-v)/2)*fac_out((J2+m2-v)/2)*fac_out((J3-J1-m2+v)/2)*fac_out((J3-J2+m1+v)/2);
		if((v/2)%2!=0)
		{
			delta=-delta;
		}
		sum=sum+1.0/delta;
	}
	delta=(J3+1)*fac_out((J1+J2-J3)/2)*fac_out((J2+J3-J1)/2)*fac_out((J3+J1-J2)/2)/fac_out((J1+J2+J3+2)/2)*fac_out((J1+m1)/2)*fac_out((J1-m1)/2)*fac_out((J2+m2)/2)*fac_out((J2-m2)/2)*fac_out((J3+m3)/2)*fac_out((J3-m3)/2);
	sum*=sqrt(delta);
	return sum;
}
inline double threej_cal(int *&j)
{
	int j1=j[0];
	int j2=j[1];
	int j3=j[2];
	int m1=j[3];
	int m2=j[4];
	int m3=j[5];
	double res=cg_cal(j1,m1,j2,m2,j3,-m3);
	if((j1-j2-m3)%4!=0)
	{
		res*=-1;
	}
	res/=sqrt(j3+1);
	return res;
}

struct node<int*,double> *threej_root=NULL;

inline int threej_cmp(int *&T1,int *&T2)
{
	return memcmp(T1,T2,sizeof(int)*6);
}
int threej_reorg(int *&j)
{
	int sign_count=0;
	int minus_count=0;
	for(int i=3;i<6;i++)
	{
		if(j[i]<0)
		{
			minus_count++;
		}
	}
	if(minus_count==2)
	{
		j[3]=-j[3];
		j[4]=-j[4];
		j[5]=-j[5];
		sign_count++;
	}
	if(minus_count==1)
	{
		int minus_i;
		int plus_i;
		for(int i=3;i<6;i++)
		{
			if(j[i]<0)
			{
				minus_i=i;
			}
			if(j[i]>0)
			{
				plus_i=i;
			}
		}
		if(j[minus_i-3]+j[minus_i]>j[plus_i-3]-j[plus_i])
		{
			j[3]=-j[3];
			j[4]=-j[4];
			j[5]=-j[5];
			sign_count++;
		}
	}
	if(j[0]>j[1]||j[0]==j[1]&&j[3]>j[4])
	{
		swap(j,j+1);
		swap(j+3,j+4);
		sign_count++;
	}
	if(j[1]>j[2]||j[1]==j[2]&&j[4]>j[5])
	{
		swap(j+1,j+2);
		swap(j+4,j+5);
		sign_count++;
	}
	if(j[0]>j[1]||j[0]==j[1]&&j[3]>j[4])
	{
		swap(j,j+1);
		swap(j+3,j+4);
		sign_count++;
	}
	if(sign_count%2==0||(j[0]+j[1]+j[2])%4==0)
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

inline double threej_out(int *j)
//the j array is deleted or stored in this function. don't try to delete it
{
	if(threej_zero(j))
	{
		delete [] j;
		return 0;
	}
	int sign=threej_reorg(j);
	return sign*tree_search_insert(threej_root,j,threej_cmp,threej_cal);
}
inline double cg_out(int j1,int m1,int j2,int m2,int j3,int m3)
{
	int *j=new int [6];
	j[0]=j1;
	j[1]=j2;
	j[2]=j3;
	j[3]=m1;
	j[4]=m2;
	j[5]=-m3;
	double res=threej_out(j);
	res*=sqrt(j3+1);
	if((j1-j2+m3)%4!=0)
	{
		res=-res;
	}
	return res;
}
struct node<int*,double> *sixj_root=NULL;

inline int sixj_cmp(int *&T1,int *&T2)
{
	return memcmp(T1,T2,sizeof(int)*6);
}



bool sixj_zero_check(int *j)
{
	if(j[2]<abs(j[0]-j[1])||j[2]>j[0]+j[1])
	{
		return(true);
	}
	if(j[0]<abs(j[4]-j[5])||j[0]>j[4]+j[5])
	{
		return(true);
	}
	if(j[1]<abs(j[3]-j[5])||j[1]>j[3]+j[5])
	{
		return(true);
	}
	if(j[2]<abs(j[3]-j[4])||j[2]>j[3]+j[4])
	{
		return(true);
	}
	return(false);
}

double delta_R(int a, int b, int c)
{
	double result;
	//int tmp1, tmp2, tmp3, tmp4;
	result=sqrt(fac_out((a+b-c)/2)*fac_out((a-b+c)/2)/fac_out((a+b+c)/2+1)*fac_out((b+c-a)/2));
	return result;
}

double sixj_cal(int *&j)
{
	double result;
	double fact;
	int tmp1=j[0]+j[1]+j[2];
	int tmp2=j[4]+j[3]+j[2];
	int tmp3=j[0]+j[4]+j[5];
	int tmp4=j[1]+j[3]+j[5];
	int k, k_min, k_max;
	int tmp5=j[0]+j[1]+j[4]+j[3];
	int tmp6=j[0]+j[3]+j[2]+j[5];
	int tmp7=j[1]+j[4]+j[2]+j[5];
	fact=delta_R(j[0], j[1], j[2])*delta_R(j[4], j[3], j[2])*delta_R(j[0], j[4], j[5])*delta_R(j[1], j[3], j[5]);
	k_min=max(tmp1, max(tmp2, max(tmp3, tmp4)));
	k_max=min(tmp5, min(tmp6, tmp7));
	result=0;
	for(k=k_min; k<=k_max; k=k+2)
	{
		if(k%4==0)
		{
			result=result+fac_out(k/2+1)/fac_out((k-tmp1)/2)/fac_out((k-tmp2)/2)/fac_out((k-tmp3)/2)/fac_out((k-tmp4)/2)/fac_out((tmp5-k)/2)/fac_out((tmp6-k)/2)/fac_out((tmp7-k)/2);
		}
		else
		{
			result=result-fac_out(k/2+1)/fac_out((k-tmp1)/2)/fac_out((k-tmp2)/2)/fac_out((k-tmp3)/2)/fac_out((k-tmp4)/2)/fac_out((tmp5-k)/2)/fac_out((tmp6-k)/2)/fac_out((tmp7-k)/2);
		}
	}
	result=result*fact;
	return result;
}

void sixj_organize(int *j)
{
	int over_case_num=0;
	if(j[0]>j[3])
	{
		over_case_num++;
	}
	if(j[1]>j[4])
	{
		over_case_num++;
	}
	if(j[2]>j[5])
	{
		over_case_num++;
	}
	if(over_case_num==1)
	{
		if(j[0]<j[3])
		{
			swap(j,j+3);
		}
		if(j[1]<j[4])
		{
			swap(j+1,j+4);
		}
		if(j[2]<j[5])
		{
			swap(j+2,j+5);
		}
	}
	if(over_case_num==2)
	{
		if(j[0]>j[3])
		{
			swap(j,j+3);
		}
		if(j[1]>j[4])
		{
			swap(j+1,j+4);
		}
		if(j[2]>j[5])
		{
			swap(j+2,j+5);
		}
	}
	if(j[2]<j[1]||((j[2]==j[1])&&(j[5]<j[4])))
	{
		swap(j+1,j+2);
		swap(j+4,j+5);
	}
	if(j[1]<j[0]||((j[1]==j[0])&&(j[4]<j[3])))
	{
		swap(j,j+1);
		swap(j+3,j+4);
	}
	if(j[2]<j[1]||((j[2]==j[1])&&(j[5]<j[4])))
	{
		swap(j+1,j+2);
		swap(j+4,j+5);
	}
}


double sixj_out(int *j)
//the j array is deleted or stored in this function. don't try to delete it
{
	if(sixj_zero_check(j))
	{
		delete [] j;
		return 0;
	}
	sixj_organize(j);
	return tree_search_insert(sixj_root,j,sixj_cmp,sixj_cal);
}	


