struct orbit
{
	int j;
	int l;
	int N;
	bool operator==(const orbit & rhs)
    {
        return( j == rhs.j) && (l == rhs.l)&&(N==rhs.N);
 
    }
};
int N_n;
int N_p;
int orbit_num_n;
int orbit_num_p;
struct orbit *orbit_n;
struct orbit *orbit_p;
size_t tb_num;
int tb_file_num;
double *ob_n;
double *ob_p;
int **tb_index;
double tb_value[2048];
int index2j(int i)
{
	if(i<orbit_num_n)
	{
		return orbit_n[i].j;
	}
	else
	{
		return orbit_p[i-orbit_num_n].j;
	}
}

struct tb
{
	int index[5];
	double value;
};

int compar(const void *p1, const void *p2)
{
	return memcmp ((void *)*((int **)p1), (void *)*((int **)p2), sizeof(int)*5);
}

void sps_int_read()
{
	tb_index=new int *[2048];

    FILE *fp;
	double temp_f;
	char filename[80]="sps_int.dat";
	if((fp=fopen(filename,"r"))==NULL)
	{
		cout<<"open file error"<<endl;
	}
	else
	{
		fscanf(fp,"%d",&N_n);
		fscanf(fp,"%d",&orbit_num_n);
		orbit_n=new struct orbit[orbit_num_n];
		ob_n=new double [orbit_num_n];
		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%d",&(orbit_n[i].j));
		}
		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%d",&(orbit_n[i].l));
		}
		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%d",&(orbit_n[i].N));
		}

		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%lf",&(ob_n[i]));
		}
		fscanf(fp,"%d",&N_p);
		fscanf(fp,"%d",&orbit_num_p);
		orbit_p=new struct orbit[orbit_num_p];
		ob_p=new double [orbit_num_p];
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%d",&(orbit_p[i].j));
		}
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%d",&(orbit_p[i].l));
		}
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%d",&(orbit_p[i].N));
		}
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%lf",&(ob_p[i]));
		}
		tb_num=0;
		for (int i = 0; i < orbit_num_n; i++)
		{
			for (int j = i; j < orbit_num_n; j++)
			{
				for (int J = abs(orbit_n[i].j - orbit_n[j].j); J <= orbit_n[i].j + orbit_n[j].j; J += 2)
				{
					if (i == j && J % 4 != 0)
					{
						continue;
					}
					for (int k = i; k < orbit_num_n; k++)
					{
						for (int l = k; l < orbit_num_n; l++)
						{
							if (i == k && j > l)
							{
								continue;
							}
							if (k == l && J % 4 != 0)
							{
								continue;
							}
							if (J < abs(orbit_n[k].j - orbit_n[l].j) || J > orbit_n[k].j + orbit_n[l].j)
							{
								continue;
							}
							if ((orbit_n[i].l + orbit_n[j].l + orbit_n[k].l + orbit_n[l].l) % 2 != 0)
							{
								continue;
							}
							tb_index[tb_num]=new int [5];
							tb_index[tb_num][0]=i+1;
							tb_index[tb_num][1]=j+1;
							tb_index[tb_num][2]=k+1;
							tb_index[tb_num][3]=l+1;
							tb_index[tb_num][4]=J/2;
							tb_num++;
						}
					}
				}
			}
		}
		for (int i = 0; i < orbit_num_p; i++)
		{
			for (int j = i; j < orbit_num_p; j++)
			{
				for (int J = abs(orbit_p[i].j - orbit_p[j].j); J <= orbit_p[i].j + orbit_p[j].j; J += 2)
				{
					if (i == j && J % 4 != 0)
					{
						continue;
					}
					for (int k = i; k < orbit_num_p; k++)
					{
						for (int l = k; l < orbit_num_p; l++)
						{
							if (i == k && j > l)
							{
								continue;
							}
							if (k == l && J % 4 != 0)
							{
								continue;
							}
							if (J < abs(orbit_p[k].j - orbit_p[l].j) || J > orbit_p[k].j + orbit_p[l].j)
							{
								continue;
							}
							if ((orbit_p[i].l + orbit_p[j].l + orbit_p[k].l + orbit_p[l].l) % 2 != 0)
							{
								continue;
							}
							tb_index[tb_num]=new int [5];
							tb_index[tb_num][0]=i+1+orbit_num_n;
							tb_index[tb_num][1]=j+1+orbit_num_n;
							tb_index[tb_num][2]=k+1+orbit_num_n;
							tb_index[tb_num][3]=l+1+orbit_num_n;
							tb_index[tb_num][4]=J/2;
							tb_num++;
						}
					}
				}
			}
		}
		for (int i = 0; i < orbit_num_n; i++)
		{
			for (int j = 0; j < orbit_num_p; j++)
			{
				for (int J = abs(orbit_n[i].j - orbit_p[j].j); J <= orbit_n[i].j + orbit_p[j].j; J += 2)
				{
					for (int k = i; k < orbit_num_n; k++)
					{
						for (int l = 0; l < orbit_num_p; l++)
						{
							if (i == k && j > l)
							{
								continue;
							}
							if (J < abs(orbit_n[k].j - orbit_p[l].j) || J > orbit_n[k].j + orbit_p[l].j)
							{
								continue;
							}
							if ((orbit_n[i].l + orbit_p[j].l + orbit_n[k].l + orbit_p[l].l) % 2 != 0)
							{
								continue;
							}
							tb_index[tb_num]=new int [5];
							tb_index[tb_num][0]=i+1;
							tb_index[tb_num][1]=j+1+orbit_num_n;
							tb_index[tb_num][2]=k+1;
							tb_index[tb_num][3]=l+1+orbit_num_n;
							tb_index[tb_num][4]=J/2;
							tb_num++;
						}
					}
				}
			}
		}
		qsort(tb_index,tb_num,sizeof(int *),compar);
		memset(tb_value,0,sizeof(double)*tb_num);
		// for(int i=0;i<tb_num;i++)
		// {
		// 	//if(fabs(tb_value[i])>1e-9)
		// 	{
		// 		cout<<tb_index[i][0]<<' '<<tb_index[i][1]<<' '<<tb_index[i][2]<<' '<<tb_index[i][3]<<' '<<tb_index[i][4]<<' '<<tb_value[i]<<endl;
		// 	}
		// }
		fscanf(fp,"%d",&tb_file_num);
		int tb_target[5];
		size_t target_index;
		int sign,j1,j2,J;
		for(int i=0;i<tb_file_num;i++)
		{
			for(int k=0;k<5;k++)
			{
				fscanf(fp,"%d",&(tb_target[k]));				
			}
			sign=1;
			J=tb_target[4]*2;
			if(tb_target[0]>tb_target[1])
			{
				swap(tb_target+0,tb_target+1);
				if(tb_target[0]>orbit_num_n)
				{
					j1=tb_target[0]-orbit_num_n-1;
					j1=orbit_p[j1].j;
				}
				else
				{
					j1=tb_target[0]-1;
					j1=orbit_n[j1].j;
				}
				if(tb_target[1]>orbit_num_n)
				{
					j2=tb_target[1]-orbit_num_n-1;
					j2=orbit_p[j2].j;
				}
				else
				{
					j2=tb_target[1]-1;
					j2=orbit_n[j2].j;
				}
				if((J-j1-j2)%4==0)
				{
					sign*=-1;
				}
			}
			if(tb_target[2]>tb_target[3])
			{
				swap(tb_target+2,tb_target+3);
				if(tb_target[2]>orbit_num_n)
				{
					j1=tb_target[2]-orbit_num_n-1;
					j1=orbit_p[j1].j;
				}
				else
				{
					j1=tb_target[2]-1;
					j1=orbit_n[j1].j;
				}
				if(tb_target[3]>orbit_num_n)
				{
					j2=tb_target[3]-orbit_num_n-1;
					j2=orbit_p[j2].j;
				}
				else
				{
					j2=tb_target[3]-1;
					j2=orbit_n[j2].j;
				}
				if((J-j1-j2)%4==0)
				{
					sign*=-1;
				}
			}
			if(tb_target[0]>tb_target[2]||(tb_target[0]==tb_target[2]&&tb_target[1]>tb_target[3]))
			{
				swap(tb_target+0,tb_target+2);
				swap(tb_target+1,tb_target+3);
			}
			if(array_search(tb_num,(void **)tb_index,(void *)tb_target,&target_index,compar))
			{
				fscanf(fp,"%lf",&temp_f);
				temp_f*=sign;
				tb_value[target_index]=temp_f;
			}
			else
			{
				fscanf(fp,"%lf",&temp_f);
			}
			
		}
        fclose(fp);
		// for(int i=0;i<tb_num;i++)
		// {
		// 	if(fabs(tb_value[i])>1e-9)
		// 	{
		// 	//	cout<<tb_index[i][0]<<' '<<tb_index[i][1]<<' '<<tb_index[i][2]<<' '<<tb_index[i][3]<<' '<<tb_index[i][4]<<' '<<tb_value[i]<<endl;
		// 	}
		// }
		// cout<<endl;
		// cin>>wat;
    }
}
int PP_J_n[2048];
int PP_pi_n[2048];
double PP_par_n[2048];
double ***PP_y_n;
int PP_num_n=0;
int QQ_num_n=0;
int QQ_kappa_n[2048];
int QQ_pi_n[2048];
int QQ_her_n[2048];
double QQ_par_n[2048];
double ***QQ_q_n;
int PP_J_p[2048];
int PP_pi_p[2048];
double PP_par_p[2048];
double ***PP_y_p;
int PP_num_p=0;
int QQ_num_p=0;
int QQ_kappa_p[2048];
int QQ_pi_p[2048];
int QQ_her_p[2048];
double QQ_par_p[2048];
double ***QQ_q_p;
int QQ_kappa_np[2048];
int QQ_pi_np[2048];
int QQ_her_np[2048];
double QQ_par_np[2048];
double ***QQ_qn_np;
double ***QQ_qp_np;
int QQ_num_np=0;

void sps_int_coll_write()
{
	FILE *fp;
	double temp_f;
	char filename[80]="sps_int_coll.dat";
	if((fp=fopen(filename,"w"))==NULL)
	{
		cout<<"open file error"<<endl;
	}
	else
	{
		fprintf(fp,"%d\n",N_n);
		fprintf(fp,"%d\n",orbit_num_n);
		for(int i=0;i<orbit_num_n;i++)
		{
			fprintf(fp,"%d ",(orbit_n[i].j));
		}
		fprintf(fp,"\n");
		for(int i=0;i<orbit_num_n;i++)
		{
			fprintf(fp,"%d ",(orbit_n[i].l));
		}
		fprintf(fp,"\n");
		for(int i=0;i<orbit_num_n;i++)
		{
			fprintf(fp,"%d ",(orbit_n[i].N));
		}
		fprintf(fp,"\n");
		for(int i=0;i<orbit_num_n;i++)
		{
			fprintf(fp,"%15.16f ",(ob_n[i]));
		}
		fprintf(fp,"\n");
		fprintf(fp,"%d\n",N_p);
		fprintf(fp,"%d\n",orbit_num_p);
		for(int i=0;i<orbit_num_p;i++)
		{
			fprintf(fp,"%d ",(orbit_p[i].j));
		}
		fprintf(fp,"\n");
		for(int i=0;i<orbit_num_p;i++)
		{
			fprintf(fp,"%d ",(orbit_p[i].l));
		}
		fprintf(fp,"\n");
		for(int i=0;i<orbit_num_p;i++)
		{
			fprintf(fp,"%d ",(orbit_p[i].N));
		}
		fprintf(fp,"\n");
		for(int i=0;i<orbit_num_p;i++)
		{
			fprintf(fp,"%15.16f ",(ob_p[i]));
		}
		fprintf(fp,"\n");
		fprintf(fp,"%d\n",PP_num_n);
		for(int i=0;i<PP_num_n;i++)
		{
			fprintf(fp,"%15.16f %d %d\n",PP_par_n[i],PP_J_n[i],PP_pi_n[i]);
			for(int p=0;p<orbit_num_n;p++)
			{
				for(int q=0;q<orbit_num_n;q++)
				{
					fprintf(fp,"%15.16f ",PP_y_n[i][p][q]);
				}
				fprintf(fp,"\n");
			}
		}
		fprintf(fp,"%d\n",QQ_num_n);
		for(int i=0;i<QQ_num_n;i++)
		{
			fprintf(fp,"%15.16f %d %d\n",QQ_par_n[i],QQ_kappa_n[i],QQ_pi_n[i]);
			for(int p=0;p<orbit_num_n;p++)
			{
				for(int q=0;q<orbit_num_n;q++)
				{
					fprintf(fp,"%15.16f ",QQ_q_n[i][p][q]);
				}
				fprintf(fp,"\n");
			}
		}
		fprintf(fp,"%d\n",PP_num_p);
		for(int i=0;i<PP_num_p;i++)
		{
			fprintf(fp,"%15.16f %d %d\n",PP_par_p[i],PP_J_p[i],PP_pi_p[i]);
			for(int p=0;p<orbit_num_p;p++)
			{
				for(int q=0;q<orbit_num_p;q++)
				{
					fprintf(fp,"%15.16f ",PP_y_p[i][p][q]);
				}
				fprintf(fp,"\n");
			}
		}
		fprintf(fp,"%d\n",QQ_num_p);
		for(int i=0;i<QQ_num_p;i++)
		{
			fprintf(fp,"%15.16f %d %d\n",QQ_par_p[i],QQ_kappa_p[i],QQ_pi_p[i]);
			for(int p=0;p<orbit_num_p;p++)
			{
				for(int q=0;q<orbit_num_p;q++)
				{
					fprintf(fp,"%15.16f ",QQ_q_p[i][p][q]);
				}
				fprintf(fp,"\n");
			}
		}

		fprintf(fp,"%d\n",QQ_num_np);
		for(int i=0;i<QQ_num_np;i++)
		{
			fprintf(fp,"%15.16f %d %d %d\n",QQ_par_np[i],QQ_kappa_np[i],QQ_pi_np[i],QQ_her_np[i]);
			for(int p=0;p<orbit_num_n;p++)
			{
				for(int q=0;q<orbit_num_n;q++)
				{
					fprintf(fp,"%15.16f ",QQ_qn_np[i][p][q]);
				}
				fprintf(fp,"\n");
			}
			for(int p=0;p<orbit_num_p;p++)
			{
				for(int q=0;q<orbit_num_p;q++)
				{
					fprintf(fp,"%15.16f ",QQ_qp_np[i][p][q]);
				}
				fprintf(fp,"\n");
			}
		}
		fclose(fp);
	}
}
void sps_int_coll_read()
{
	FILE *fp;
	double temp_f;
	char filename[80]="sps_int_coll.dat";
	if((fp=fopen(filename,"r"))==NULL)
	{
		cout<<"open file error"<<endl;
	}
	else
	{
		fscanf(fp,"%d",&N_n);
		fscanf(fp,"%d",&orbit_num_n);
		orbit_n=new struct orbit[orbit_num_n];
		ob_n=new double [orbit_num_n];
		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%d",&(orbit_n[i].j));
		}
		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%d",&(orbit_n[i].l));
		}
		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%d",&(orbit_n[i].N));
		}

		for(int i=0;i<orbit_num_n;i++)
		{
			fscanf(fp,"%lf",&(ob_n[i]));
		}
		fscanf(fp,"%d",&N_p);
		fscanf(fp,"%d",&orbit_num_p);
		orbit_p=new struct orbit[orbit_num_p];
		ob_p=new double [orbit_num_p];
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%d",&(orbit_p[i].j));
		}
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%d",&(orbit_p[i].l));
		}
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%d",&(orbit_p[i].N));
		}
		for(int i=0;i<orbit_num_p;i++)
		{
			fscanf(fp,"%lf",&(ob_p[i]));
		}
		fscanf(fp,"%d",&PP_num_n);
		PP_y_n=new double **[PP_num_n];
		for(int i=0;i<PP_num_n;i++)
		{
			PP_y_n[i]=new double *[orbit_num_n];
			for(int p=0;p<orbit_num_n;p++)
			{
				PP_y_n[i][p]=new double [orbit_num_n];
			}
		}
		for(int i=0;i<PP_num_n;i++)
		{
			fscanf(fp,"%lf %d %d",PP_par_n+i,PP_J_n+i,PP_pi_n+i);
			for(int p=0;p<orbit_num_n;p++)
			{
				for(int q=0;q<orbit_num_n;q++)
				{
					fscanf(fp,"%lf ",PP_y_n[i][p]+q);
				}
			}
		}
		fscanf(fp,"%d",&QQ_num_n);
		QQ_q_n=new double **[QQ_num_n];
		for(int i=0;i<QQ_num_n;i++)
		{
			QQ_q_n[i]=new double *[orbit_num_n];
			for(int p=0;p<orbit_num_n;p++)
			{
				QQ_q_n[i][p]=new double [orbit_num_n];
			}
		}
		for(int i=0;i<QQ_num_n;i++)
		{
			fscanf(fp,"%lf %d %d %d\n",QQ_par_n+i,QQ_kappa_n+i,QQ_pi_n+i,QQ_her_n+i);
			for(int p=0;p<orbit_num_n;p++)
			{
				for(int q=0;q<orbit_num_n;q++)
				{
					fscanf(fp,"%lf ",QQ_q_n[i][p]+q);
				}
			}
		}
		fscanf(fp,"%d",&PP_num_p);
		PP_y_p=new double **[PP_num_p];
		for(int i=0;i<PP_num_p;i++)
		{
			PP_y_p[i]=new double *[orbit_num_p];
			for(int p=0;p<orbit_num_p;p++)
			{
				PP_y_p[i][p]=new double [orbit_num_p];
			}
		}
		for(int i=0;i<PP_num_p;i++)
		{
			fscanf(fp,"%lf %d %d",PP_par_p+i,PP_J_p+i,PP_pi_p+i);
			for(int p=0;p<orbit_num_p;p++)
			{
				for(int q=0;q<orbit_num_p;q++)
				{
					fscanf(fp,"%lf ",PP_y_p[i][p]+q);
				}
			}
		}
		fscanf(fp,"%d",&QQ_num_p);
		QQ_q_p=new double **[QQ_num_p];
		for(int i=0;i<QQ_num_p;i++)
		{
			QQ_q_p[i]=new double *[orbit_num_p];
			for(int p=0;p<orbit_num_p;p++)
			{
				QQ_q_p[i][p]=new double [orbit_num_p];
			}
		}
		for(int i=0;i<QQ_num_p;i++)
		{
			fscanf(fp,"%lf %d %d %d\n",QQ_par_p+i,QQ_kappa_p+i,QQ_pi_p+i,QQ_her_p+i);
			for(int p=0;p<orbit_num_p;p++)
			{
				for(int q=0;q<orbit_num_p;q++)
				{
					fscanf(fp,"%lf ",QQ_q_p[i][p]+q);
				}
			}
		}
		fscanf(fp,"%d",&QQ_num_np);
		QQ_qn_np=new double **[QQ_num_np];
		QQ_qp_np=new double **[QQ_num_np];
		for(int i=0;i<QQ_num_np;i++)
		{
			QQ_qn_np[i]=new double *[orbit_num_n];
			for(int p=0;p<orbit_num_n;p++)
			{
				QQ_qn_np[i][p]=new double [orbit_num_n];
			}
			QQ_qp_np[i]=new double *[orbit_num_p];
			for(int p=0;p<orbit_num_p;p++)
			{
				QQ_qp_np[i][p]=new double [orbit_num_p];
			}
		}
		for(int i=0;i<QQ_num_np;i++)
		{
			fscanf(fp,"%lf %d %d %d\n",QQ_par_np+i,QQ_kappa_np+i,QQ_pi_np+i,QQ_her_np+i);
			for(int p=0;p<orbit_num_n;p++)
			{
				for(int q=0;q<orbit_num_n;q++)
				{
					fscanf(fp,"%lf ",QQ_qn_np[i][p]+q);
				}
			}
			for(int p=0;p<orbit_num_p;p++)
			{
				for(int q=0;q<orbit_num_p;q++)
				{
					fscanf(fp,"%lf ",QQ_qp_np[i][p]+q);
				}
			}
		}
		fclose(fp);
	}
}