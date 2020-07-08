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
clock_t t1,t2;
cycles_t t11,t22,t3,t4;
cycles_t t_sum=0;
int worr;
char yorn;
char yorn_eig_display[5];
int M_tot;
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
//#include "ham_q_mat_out.cpp"
#include "basis_tot_out.cpp"
#include "vec_read.cpp"
#include "tran_input_read.cpp"
#include "tran_q_mat_out.cpp"
#include "tran_cal.cpp"
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
	
	cout<<"Q_M mat index start from "<<INT_MAX/2<<endl;
	cout<<"sp vector index start from "<<INT_MAX/4<<endl;
	
	sps_read();
	PQ0_add_init();
	J_pm_out(0);
	J_pm_out(1);
	P_read();
	vec_read();
	basis_out(0);
	basis_out(1);
	block_mat_init();
	ham_read();
	tran_input_read();
	pow_gamma_out(0);
	pow_q_out(0);
	pow_gamma_out(1);
	pow_q_out(1);

	block_mat_read();
	basis_tot_out(M_tot);
	
	cout<<"are transition Q matrix files avilable?(y/n):"<<endl;
	cout<<"files avilable only under modification of interesting tran_mat elements (modification of interesting tran Q operator is not the case), or the M_tot modification in previous spectra.out operation"<<endl;
	cin>>yorn;
	//stop here**** q_mat_write_read system*********
	tran_q_block_init();
	if(yorn=='y')
	{
		tran_q_block_read();
	}
	t1=clock();
	ove_tail_init();
	for(int norp=0;norp<2;norp++)
	{
		for(int pi_l=0;pi_l<2;pi_l++)
		{
			for(int pi_r=0;pi_r<2;pi_r++)
			{
				tran_q_block_out(norp,pi_l,pi_r);
			}
		}
	}
	t2=clock();
	tran_q_block_write();
	double t_mat=double(t2-t1)/CLOCKS_PER_SEC;
	t1=clock();
	tran_cal();
	t2=clock();
	double t_tran=double(t2-t1)/CLOCKS_PER_SEC;
	cout<<"t_mat="<<t_mat<<" t_tran="<<t_tran<<endl;
	return 0;
}