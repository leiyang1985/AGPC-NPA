#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<signal.h>
#include <string.h>
#include<vector>
#include<map>
#include <algorithm> 
#include <utility>
#include<time.h>
#include <stdio.h> /* printf */ 
#include <unistd.h> /* alarm, pause */ 
#include <sys/types.h> 
#include <signal.h> /* signal,kill */
#include<linux/types.h>
using namespace std;
#include "mkl.h"
//以上都是系统或标准库中的引用，你的不用管



#include "lib.h"//我自已写的一些常用函数，你的不用关心
#include "am.cpp"//我自已写的角动量代数相关的函数，主要是程序中要用到fac_out()函数，是阶乘
#include "ove_for.cpp"//算overlap的函数文件



int main()
{
	int N=3;
	cin>>N;//N为对的数目
	ove_for_init(N);//对N作划分，划分结果放入nt_store中
	// cout<<"rs="<<endl;
	// for(int i=0;i<rs_store.size();i++)
	// {
	// 	for(int k=0;k<N;k++)
	// 	{
	// 		cout<<int(rs_store[i][k])<<' ';
	// 	}
	// 	cout<<endl;
	// }
	// cout<<"nt="<<endl;
	// for(int i=0;i<nt_store.size();i++)
	// {
	// 	for(int k=0;k<nt_store[i].size();k++)
	// 	{
	// 		cout<<int(nt_store[i][k])<<' ';
	// 	}
	// 	cout<<endl;
	// }
	//cout<<"r="<<endl;
	// clock_t t,t1,t2;
	// t1=clock();
	
	// cout<<ove_for_cal(N)<<endl;//开始算
	// t2=clock();
	// t=t2-t1;
	// cout<<double(t)/CLOCKS_PER_SEC<<endl;
}
