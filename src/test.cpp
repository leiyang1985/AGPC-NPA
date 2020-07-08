#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<signal.h>
#include<vector>
#include<array>
#include <algorithm> 
#include <utility>
#include <string.h>
#include<limits.h>
#include<stddef.h>
#include<time.h>
#include <unistd.h> /* alarm, pause */ 
#include <sys/types.h> 
#include <signal.h> /* signal,kill */
#include<linux/types.h>
using namespace std;
int wat;
bool flag_print=false;
size_t found=0;
size_t all=0;


#include "mkl.h"
#include "lib.h"//for linear operation
clock_t t,t1,t2;
cycles_t t11,t22,t3,t4;
cycles_t t_sum=0;
int worr;
int M_tot=0;
char yorn;
char yorn_eig_display[5];
FILE *fp;
//#include "array_operation.cpp"//memory operation
#include "am.cpp"
//#include "struct.cpp"
#include "sps_read.cpp"
#include "PQ_operation.cpp"
#include "P_read.cpp"
#include "basis_out.cpp"
#include "mm.cpp"
#include "tr.cpp"
#include "PP_dot.cpp"
#include "CC_dot.cpp"
#include "CC_con.cpp"
#include "PP_con.cpp"
#include "PQ_con.cpp"
#include "CQ_con.cpp"

#include "overlap.cpp"
//#include "P_npa_read.cpp"
#include "Qmu_out.cpp"
#include "qq.cpp"
#include "q.cpp"
#include "pp.cpp"
#include "ham_read.cpp"
#include "block_mat_out.cpp"
#include "ham_q_mat_out.cpp"
#include "basis_tot_out.cpp"
#include "eig.cpp"
int main()
{
	// cout<<cg_cal(11,-9,11,7,4,-2)<<endl;
	// cin>>wat;
	mkl_set_num_threads(1);	
	// int *j=new int [6];
	// while(true)
	// {
	// 	j=new int [6];
	// 	for(int i=0;i<6;i++)
	// 	{
	// 		cin>>j[i];
	// 	}
	// 	cout<<sixj_out(j)<<endl;
	// }
	t1=clock();
	cout<<"Q_M mat index start from "<<INT_MAX/2<<endl;
	cout<<"sp vector index start from "<<INT_MAX/4<<endl;
	
	sps_read();
	PQ0_add_init();
	J_pm_out(0);
	J_pm_out(1);
	P_read();
	cout<<"input the M_tot(M_tot=0 disables some transition calculation)=";
	cin>>M_tot;
	if((M_tot+N[0]+N[1])%2!=0)
	{
		cout<<"illegal of inputted M_tot for the particle number specified in sps.dat. Abort"<<endl;
		return 0;
	}
	basis_out(0);
	basis_out(1);
	block_mat_init();
	ham_read();
	pow_gamma_out(0);
	pow_q_out(0);
	pow_gamma_out(1);
	pow_q_out(1);
	cout<<"are previous matrix files avilable?(y/n):"<<endl;
	cout<<"files avilable only under the modification of M_tot or Ham parameters."<<endl;
	cout<<"by inputting n, the previous binary files will be erased no matter."<<endl;
	cin>>yorn;
	if(yorn=='y')
	{
		block_mat_read();
	}
	else
	{
		system("rm *.bin");
	}
	
	cout<<"please specify which expectaions are displayed in sequence of"<<endl;
	cout<<"pair(unpaired paritcl) num, sp occupation, PP, QQ_norp, QQ_np (y/n):"<<endl;
	
	for(int k=0;k<5;k++)
	{
		cin>>yorn_eig_display[k];
	}
	ove_tail_init();

	t1=clock();
	for(int norp=0;norp<2;norp++)
	{
		for(int pi=0;pi<2;pi++)
		{
			basis_filter(norp,pi);
			block_mat_out(norp,pi);
		}
	}	
	block_mat_write();
	ham_q_block_init();
	if(yorn=='y')
	{
		ham_q_block_read();
	}
	for(int norp=0;norp<2;norp++)
	{
		for(int pi_l=0;pi_l<2;pi_l++)
		{
			for(int pi_r=0;pi_r<2;pi_r++)
			{
				ham_q_block_out(norp,pi_l,pi_r);
			}
		}
	}
	ham_q_block_write();
	t2=clock();
	t=t2-t1;
	double t_mat_cal=double(t)/CLOCKS_PER_SEC;
	t1=clock();
	basis_tot_out(M_tot);
	eig_input_read();
	dig();
	t2=clock();
	t=t2-t1;
	double t_dig=double(t)/CLOCKS_PER_SEC;
	cout<<"t for mat cal="<<t_mat_cal<<"s, for dig="<<t_dig<<"s"<<endl;
	return 0;
}