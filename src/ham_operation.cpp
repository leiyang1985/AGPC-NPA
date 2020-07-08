int Fk_num_n;
int **Fk_index_n;
double Fk_value_n[2048];
void Fk_index_evalue_n()
{
	Fk_num_n = 0;
	Fk_index_n = new int *[2048];
	for (int i = 0; i < orbit_num_n; i++)
	{
		for (int j = i; j < orbit_num_n; j++)
		{
			for (int J = abs(orbit_n[i].j - orbit_n[j].j); J <= orbit_n[i].j + orbit_n[j].j; J += 2)
			{
				for (int k = i; k < orbit_num_n; k++)
				{
					for (int l = k; l < orbit_num_n; l++)
					{
						if (i == k && j > l)
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
						Fk_index_n[Fk_num_n] = new int[5];
						Fk_index_n[Fk_num_n][0] = i + 1;
						Fk_index_n[Fk_num_n][1] = j + 1;
						Fk_index_n[Fk_num_n][2] = k + 1;
						Fk_index_n[Fk_num_n][3] = l + 1;
						Fk_index_n[Fk_num_n][4] = J / 2;
						Fk_num_n++;
					}
				}
			}
		}
	}
	qsort(Fk_index_n, Fk_num_n, sizeof(int *), compar);
	memset(Fk_value_n, 0, sizeof(double) * Fk_num_n);
}

int tb_num_n;
int **tb_index_n;
double tb_value_n[2048];
double tb_P_value_n[2048];
double tb_Q_value_n[2048];
void tb_evalue_n() //it also inintal tb_value
{
	tb_num_n = 0;
	tb_index_n = new int *[2048];
	for (int i = 0; i < tb_num; i++)
	{
		if (tb_index[i][0] < orbit_num_n + 1)
		{
			if (tb_index[i][1] < orbit_num_n + 1)
			{
				if (tb_index[i][2] < orbit_num_n + 1)
				{
					if (tb_index[i][3] < orbit_num_n + 1)
					{
						tb_index_n[tb_num_n] = new int[5];
						memcpy(tb_index_n[tb_num_n], tb_index[i], sizeof(int) * 5);
						tb_value_n[tb_num_n] = tb_value[i];
						tb_num_n++;
					}
				}
			}
		}
	}
}
struct PQ_eig_vec
{
	double eig;
	double *vec;
	int coll_index;
};

struct PQ_eig_vec *P_eig_vec_n;
double *P_eig_n;
int P_eig_num_n = 0;
int P_eig_vec_num_n = 0;
int *P_eig_order_n;

int compar_P_eig_vec_n(const void *p1, const void *p2)
{
	double eig1 = fabs(P_eig_vec_n[*((int *)p1)].eig);
	double eig2 = fabs(P_eig_vec_n[*((int *)p2)].eig);
	if (eig1 > eig2)
	{
		return -1;
	}
	else if (eig1 < eig2)
	{
		return 1;
	}
	return 0;
}
void dig_P_init_n()
{
	int dim;
	P_eig_vec_num_n = 0;
	for (int i = 0; i < P_coll_stru_num_n; i++)
	{
		dim = P_coll_stru_tab_n[i].dim;
		P_eig_vec_num_n += dim;
	}
	P_eig_vec_n = new struct PQ_eig_vec[P_eig_vec_num_n];
	P_eig_order_n = new int[P_eig_vec_num_n];
	P_eig_n = new double[P_eig_vec_num_n];
	P_eig_num_n = P_eig_vec_num_n;
	P_eig_vec_num_n = 0;
	for (int i = 0; i < P_coll_stru_num_n; i++)
	{
		dim = P_coll_stru_tab_n[i].dim;
		for (int k = 0; k < dim; k++)
		{
			P_eig_vec_n[P_eig_vec_num_n].vec = new double[dim];
			P_eig_vec_n[P_eig_vec_num_n].coll_index = i;
			P_eig_order_n[P_eig_vec_num_n] = P_eig_vec_num_n;
			P_eig_vec_num_n++;
		}
	}
}
void dig_P_eig_vec_n()
{
	int target[5];
	size_t target_index;
	char a = 'V';
	char b = 'U';
	int info = 0;
	double *h;
	double *eig;
	double *work;
	int lwork;
	int dim;
	int dim2;
	int inc = 1;
	P_eig_vec_num_n = 0;
	for (int i = 0; i < P_coll_stru_num_n; i++)
	{
		dim = P_coll_stru_tab_n[i].dim;
		dim2 = dim * dim;
		h = new double[dim2];
		eig = new double[dim];
		memset(h, 0, sizeof(double) * dim2);
		lwork = dim * 3;
		work = new double[lwork];
		target[4] = P_coll_stru_tab_n[i].J / 2;
		for (int k = 0; k < dim; k++)
		{
			for (int l = k; l < dim; l++)
			{
				target[0] = P_coll_stru_tab_n[i].i[k] + 1;
				target[1] = P_coll_stru_tab_n[i].j[k] + 1;
				target[2] = P_coll_stru_tab_n[i].i[l] + 1;
				target[3] = P_coll_stru_tab_n[i].j[l] + 1;

				if (array_search(tb_num_n, (void **)tb_index_n, (void *)target, &target_index, compar))
				{
					h[k + l * dim] += tb_P_value_n[target_index];
				}
			}
		}
		bool flag_dia_able = false;
		for (int k = 0; k < dim2; k++)
		{
			if (fabs(h[k]) > 0)
			{
				flag_dia_able = true;
				break;
			}
		}
		if (flag_dia_able)
		{
			if (flag_print)
			{
				cout << "J=" << P_coll_stru_tab_n[i].J << endl;
				for (int k = 0; k < dim; k++)
				{
					for (int l = 0; l < dim; l++)
					{
						cout << h[k + l * dim] << ' ';
					}
					cout << endl;
				}
			}
			dsyev_(&a, &b, &dim, h, &dim, eig, work, &lwork, &info);
			if (flag_print)
			{
				cout << "eig=";
				for (int k = 0; k < dim; k++)
				{
					cout << eig[k] << ' ';
				}
				cout << endl;
			}
			for (int k = 0; k < dim; k++)
			{
				P_eig_vec_n[P_eig_vec_num_n].eig = eig[k];
				memcpy(P_eig_vec_n[P_eig_vec_num_n].vec, h + k * dim, sizeof(double) * dim);
				P_eig_vec_num_n++;
			}
		}
		else
		{
			for (int k = 0; k < dim; k++)
			{
				P_eig_vec_n[P_eig_vec_num_n].eig = 0;
				P_eig_vec_n[P_eig_vec_num_n].coll_index = i;
				memset(P_eig_vec_n[P_eig_vec_num_n].vec, 0, sizeof(double) * dim);
				P_eig_vec_num_n++;
			}
		}

		delete[] h;
		delete[] eig;
		delete[] work;
	}
}
/********************below is the one actually usefull P eig******************************/

void dig_P_out_n()
{
	int target[5];
	size_t target_index;
	char a = 'V';
	char b = 'U';
	int info = 0;
	double *h;
	double *eig;
	double *work;
	int lwork;
	int dim;
	int dim2;
	int inc = 1;
	double y[orbit_num_n][orbit_num_n];
	PP_y_n=new double **[2048];
	for(int i=0;i<2048;i++)
	{
		PP_y_n[i]=new double *[orbit_num_n];
		for(int p=0;p<orbit_num_n;p++)
		{
			PP_y_n[i][p]=new double [orbit_num_n];
		}
	}
	PP_num_n=0;
	for (int i = 0; i < P_coll_stru_num_n; i++)
	{
		dim = P_coll_stru_tab_n[i].dim;
		dim2 = dim * dim;
		h = new double[dim2];
		eig = new double[dim];
		memset(h, 0, sizeof(double) * dim2);
		lwork = dim * 3;
		work = new double[lwork];
		target[4] = P_coll_stru_tab_n[i].J / 2;
		for (int k = 0; k < dim; k++)
		{
			for (int l = k; l < dim; l++)
			{
				target[0] = P_coll_stru_tab_n[i].i[k] + 1;
				target[1] = P_coll_stru_tab_n[i].j[k] + 1;
				target[2] = P_coll_stru_tab_n[i].i[l] + 1;
				target[3] = P_coll_stru_tab_n[i].j[l] + 1;

				if (array_search(tb_num, (void **)tb_index, (void *)target, &target_index, compar))
				{
					h[k + l * dim] += tb_value[target_index];
				}
			}
		}
		bool flag_dia_able = false;
		for (int k = 0; k < dim2; k++)
		{
			if (fabs(h[k]) > 0)
			{
				flag_dia_able = true;
				break;
			}
		}
		if (flag_dia_able)
		{
			cout << "J/pi=" << P_coll_stru_tab_n[i].J << ' '<<P_coll_stru_tab_n[i].pi<<endl;
			dsyev_(&a, &b, &dim, h, &dim, eig, work, &lwork, &info);
			for (int k = 0; k < dim; k++)
			{
				if(fabs(eig[k])<1e-3)
				{
					continue;
				}
				PP_par_n[PP_num_n]=eig[k];
				PP_J_n[PP_num_n]=P_coll_stru_tab_n[i].J;
				PP_pi_n[PP_num_n]=P_coll_stru_tab_n[i].pi;
				cout << "eig= "<<eig[k]<<endl;
				cout<<"vec= ";
				for(int p=0;p<dim;p++)
				{
					cout<<h[k*dim+p]<<' ';
				}
				cout<<endl;
				for(int p=0;p<orbit_num_n;p++)
				{
					for(int q=0;q<orbit_num_n;q++)
					{
						y[p][q]=0;
						
					}
				}
				for(int veci=0;veci<dim;veci++)
				{
					int ii=P_coll_stru_tab_n[i].i[veci];
					int jj=P_coll_stru_tab_n[i].j[veci];
					if(ii==jj)
					{
						y[ii][jj]+=h[k*dim+veci]/sqrt(2);
					}
					else
					{
						y[ii][jj]+=h[k*dim+veci]/2;
						if((orbit_n[ii].j+orbit_n[jj].j+P_coll_stru_tab_n[i].J)%4==0)
						{
							y[jj][ii]-=h[k*dim+veci]/2;
						}
						else
						{
							y[jj][ii]+=h[k*dim+veci]/2;
						}
						
					}
				}
				cout<<"P="<<endl;
				for(int p=0;p<orbit_num_n;p++)
				{
					for(int q=0;q<orbit_num_n;q++)
					{
						cout<<y[p][q]<<' ';
						PP_y_n[PP_num_n][p][q]=y[p][q];
					}
					cout<<endl;
				}
				cout<<endl;
				PP_num_n++;
			}
		}

		delete[] h;
		delete[] eig;
		delete[] work;
	}
}

void dig_P_out_p()
{
	int target[5];
	size_t target_index;
	char a = 'V';
	char b = 'U';
	int info = 0;
	double *h;
	double *eig;
	double *work;
	int lwork;
	int dim;
	int dim2;
	int inc = 1;
	double y[orbit_num_p][orbit_num_p];
	PP_y_p=new double **[2048];
	for(int i=0;i<2048;i++)
	{
		PP_y_p[i]=new double *[orbit_num_p];
		for(int p=0;p<orbit_num_p;p++)
		{
			PP_y_p[i][p]=new double [orbit_num_p];
		}
	}
	PP_num_p=0;
	for (int i = 0; i < P_coll_stru_num_p; i++)
	{
		dim = P_coll_stru_tab_p[i].dim;
		dim2 = dim * dim;
		h = new double[dim2];
		eig = new double[dim];
		memset(h, 0, sizeof(double) * dim2);
		lwork = dim * 3;
		work = new double[lwork];
		target[4] = P_coll_stru_tab_p[i].J / 2;
		for (int k = 0; k < dim; k++)
		{
			for (int l = k; l < dim; l++)
			{
				target[0] = P_coll_stru_tab_p[i].i[k] + 1+orbit_num_n;
				target[1] = P_coll_stru_tab_p[i].j[k] + 1+orbit_num_n;
				target[2] = P_coll_stru_tab_p[i].i[l] + 1+orbit_num_n;
				target[3] = P_coll_stru_tab_p[i].j[l] + 1+orbit_num_n;

				if (array_search(tb_num, (void **)tb_index, (void *)target, &target_index, compar))
				{
					h[k + l * dim] += tb_value[target_index];
				}
			}
		}
		bool flag_dia_able = false;
		for (int k = 0; k < dim2; k++)
		{
			if (fabs(h[k]) > 0)
			{
				flag_dia_able = true;
				break;
			}
		}
		if (flag_dia_able)
		{
			cout << "J/pi=" << P_coll_stru_tab_p[i].J << ' '<<P_coll_stru_tab_p[i].pi<<endl;
			dsyev_(&a, &b, &dim, h, &dim, eig, work, &lwork, &info);
			for (int k = 0; k < dim; k++)
			{
				if(fabs(eig[k])<1e-3)
				{
					continue;
				}
				PP_par_p[PP_num_p]=eig[k];
				PP_J_p[PP_num_p]=P_coll_stru_tab_p[i].J;
				PP_pi_p[PP_num_p]=P_coll_stru_tab_p[i].pi;
				cout << "eig= "<<eig[k]<<endl;
				cout<<"vec= ";
				for(int p=0;p<dim;p++)
				{
					cout<<h[k*dim+p]<<' ';
				}
				cout<<endl;
				for(int p=0;p<orbit_num_p;p++)
				{
					for(int q=0;q<orbit_num_p;q++)
					{
						y[p][q]=0;
						
					}
				}
				for(int veci=0;veci<dim;veci++)
				{
					int ii=P_coll_stru_tab_p[i].i[veci];
					int jj=P_coll_stru_tab_p[i].j[veci];
					if(ii==jj)
					{
						y[ii][jj]+=h[k*dim+veci]/sqrt(2);
					}
					else
					{
						y[ii][jj]+=h[k*dim+veci]/2;
						if((orbit_p[ii].j+orbit_p[jj].j+P_coll_stru_tab_p[i].J)%4==0)
						{
							y[jj][ii]-=h[k*dim+veci]/2;
						}
						else
						{
							y[jj][ii]+=h[k*dim+veci]/2;
						}
						
					}
				}
				cout<<"P="<<endl;
				for(int p=0;p<orbit_num_p;p++)
				{
					for(int q=0;q<orbit_num_p;q++)
					{
						cout<<y[p][q]<<' ';
						PP_y_p[PP_num_p][p][q]=y[p][q];
					}
					cout<<endl;
				}
				cout<<endl;
				PP_num_p++;
			}
		}

		delete[] h;
		delete[] eig;
		delete[] work;
	}
}

void svd_Q_out_np()
{
	int target[5];
	size_t target_index;
	char a = 'V';
	char b = 'U';
	int info = 0;
	double *h;
	double *eig;
	double *work;
	double *u;
	double *vt;
	int lwork;
	int dim_n;
	int dim_p;
	int dim_min;
	int dim_max;
	int dim2;
	int inc = 1;
	int kappa;
	int pi;
	int Qp_index;
	double delta;
	int orbit_num_max;
	if(orbit_num_n>orbit_num_p)
	{
		orbit_num_max=orbit_num_n;
	}
	else
	{
		orbit_num_max=orbit_num_p;
	}
	double tb_value_test[tb_num];
	memset(tb_value_test,0,sizeof(double)*tb_num);
	double y[orbit_num_max][orbit_num_max];
	QQ_num_np=0;
	QQ_qn_np=new double **[2048];
	QQ_qp_np=new double **[2048];
	for(int i=0;i<2048;i++)
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

	for (int i = 0; i < Q_coll_stru_num_n; i++)
	{
		dim_n = Q_coll_stru_tab_n[i].dim;
		kappa=Q_coll_stru_tab_n[i].J;
		pi=Q_coll_stru_tab_n[i].pi;
		int Qp_index=JP_in_PQ_coll_stru_tab(kappa,pi,Q_coll_stru_num_p,Q_coll_stru_tab_p);
		if(Qp_index==-1)
		{
			continue;
		}
		dim_p=Q_coll_stru_tab_p[Qp_index].dim;
		if(dim_n<dim_p)
		{
			dim_min=dim_n;
			dim_max=dim_p;
		}
		else
		{
			dim_min=dim_p;
			dim_max=dim_n;
		}
		h = new double[dim_n*dim_p];
		eig = new double[dim_min];
		memset(h, 0, sizeof(double) * dim_n*dim_p);
		
		for (int k = 0; k < dim_n; k++)
		{
			for (int l = 0; l < dim_p; l++)
			{
				int ii=Q_coll_stru_tab_n[i].i[k];
				int jj=Q_coll_stru_tab_n[i].j[k];
				int kk=Q_coll_stru_tab_p[Qp_index].i[l];
				int ll=Q_coll_stru_tab_n[Qp_index].j[l];
				target[0] = ii + 1;
				target[1] = jj + 1;
				target[2] = kk + 1+orbit_num_n;
				target[3] = ll + 1+orbit_num_n;
				swap(target+1,target+2);
				if(target[0]>target[2]||(target[0]==target[2]&&target[1]>target[3]))
				{
					swap(target,target+2);
					swap(target+1,target+3);
				}
				for(target[4]=abs(orbit_n[target[0]-1].j-orbit_p[target[1]-1-orbit_num_n].j)/2;target[4]<=(orbit_n[target[0]-1].j+orbit_p[target[1]-1-orbit_num_n].j)/2;target[4]++)
				{
					if(target[4]<abs(orbit_n[target[2]-1].j-orbit_p[target[3]-1-orbit_num_n].j)/2||target[4]>(orbit_n[target[2]-1].j+orbit_p[target[3]-1-orbit_num_n].j)/2)
					{
						continue;
					}
					if (array_search(tb_num, (void **)tb_index, (void *)target, &target_index, compar))
					{
						delta=tb_value[target_index];
						if((orbit_n[jj].j+orbit_p[kk].j+target[4]*2)%4!=0)
						{
							delta*=-1;
						}
						delta*=target[4]*2+1;
						delta*=sixj_out(orbit_n[ii].j,orbit_p[kk].j,target[4]*2,
						orbit_p[ll].j,orbit_n[jj].j,kappa);
						h[k + l * dim_n] += delta;
					}
				}
			}
		}
		bool flag_dia_able = false;
		for (int k = 0; k < dim_n*dim_p; k++)
		{
			if (fabs(h[k]) > 0)
			{
				flag_dia_able = true;
				break;
			}
		}
		if (flag_dia_able)
		{
			cout << "kappa/pi=" << kappa << ' '<<pi<<endl;
			int lda=dim_n;
			int lwork=-1;
			char jobu='S';
			char jobvt='S';
			u=new double [dim_n*dim_min];
			int ldu=dim_n;
			vt=new double [dim_min*dim_p];
			int ldvt=dim_min;
			work=new double [1];
			int info;
			dgesvd_(&jobu,&jobvt,&dim_n,&dim_p,h,&lda,eig,u,&ldu,vt,&ldvt,work,&lwork,&info);
			lwork=work[0];
			delete [] work;
			work=new double [lwork];
			dgesvd_(&jobu,&jobvt,&dim_n,&dim_p,h,&lda,eig,u,&ldu,vt,&ldvt,work,&lwork,&info);
			mkl_dimatcopy('c','t',dim_min,dim_p,1,vt,dim_min,dim_p);
			for(int k=0;k<dim_min;k++)
			{
				cout<<eig[k]<<' ';
			}
			cout<<endl;
			for (int k = 0; k < dim_min; k++)
			{
				if(fabs(eig[k])<1e-3)
				{
					continue;
				}
				QQ_par_np[QQ_num_np]=eig[k];
				QQ_kappa_np[QQ_num_np]=kappa;
				QQ_pi_np[QQ_num_np]=pi;
				for(int p=0;p<dim_n;p++)
				{
					for(int q=0;q<dim_p;q++)
					{
						delta=eig[k]*u[p+k*dim_n]*vt[q+k*dim_p];
						int ii=Q_coll_stru_tab_n[i].i[p];
						int jj=Q_coll_stru_tab_n[i].j[p];
						int kk=Q_coll_stru_tab_p[Qp_index].i[q];
						int ll=Q_coll_stru_tab_p[Qp_index].j[q];
						for(int J=abs(orbit_n[ii].j-orbit_p[kk].j);J<=orbit_n[ii].j+orbit_p[kk].j;J+=2)
						{
							double delta1;
							if(J<abs(orbit_n[jj].j-orbit_p[ll].j)||J>orbit_n[jj].j+orbit_p[ll].j)
							{
								continue;
							}
							if((orbit_n[jj].j+orbit_p[kk].j+J)%4!=0)
							{
								delta1=-delta;
							}
							else
							{
								delta1=delta;
							}
							delta1*=kappa+1;
							delta1*=sixj_out(orbit_n[ii].j,orbit_n[jj].j,kappa,
											orbit_p[ll].j,orbit_p[kk].j,J);
							target[0]=ii+1;
							target[1]=kk+1+orbit_num_n;
							target[2]=jj+1;
							target[3]=ll+1+orbit_num_n;
							target[4]=J/2;
							if (array_search(tb_num, (void **)tb_index, (void *)target, &target_index, compar))
							{
								tb_value_test[target_index]+=delta1;
							}
						}
					}
				}
				cout << "eig= "<<eig[k]<<endl;
				cout<<"vecl= ";
				for(int p=0;p<dim_n;p++)
				{
					cout<<u[p+k*dim_n]<<' ';
				}
				cout<<endl;
				for(int p=0;p<orbit_num_n;p++)
				{
					for(int q=0;q<orbit_num_n;q++)
					{
						y[p][q]=0;
					}
				}
				for(int veci=0;veci<dim_n;veci++)
				{
					int ii=Q_coll_stru_tab_n[i].i[veci];
					int jj=Q_coll_stru_tab_n[i].j[veci];
					y[ii][jj]+=u[veci+k*dim_n];
				}
				bool isher=true;
				bool isantiher=true;
				cout<<"Qn="<<endl;
				for(int p=0;p<orbit_num_n;p++)
				{
					for(int q=0;q<orbit_num_n;q++)
					{
						cout<<y[p][q]<<' ';
						QQ_qn_np[QQ_num_np][p][q]=y[p][q];
					}
					cout<<endl;
				}
				for(int p=0;p<orbit_num_n;p++)
				{
					for(int q=p;q<orbit_num_n;q++)
					{
						double delta=y[p][q];
						if((orbit_n[p].j-orbit_n[q].j)%4==0)
						{
							delta-=y[q][p];
						}
						else
						{
							delta+=y[q][p];
						}
						
						if(fabs(delta)>1e-7)
						{
							isher=false;
						}
						delta=y[p][q];
						if((orbit_n[p].j+orbit_n[q].j)%4==0)
						{
							delta-=y[q][p];
						}
						else
						{
							delta+=y[q][p];
						}
						
						if(fabs(delta)>1e-7)
						{
							isantiher=false;
						}
					}
				}
				cout<<"vecr= ";
				for(int p=0;p<dim_p;p++)
				{
					cout<<vt[p+k*dim_p]<<' ';
				}
				cout<<endl;
				for(int p=0;p<orbit_num_p;p++)
				{
					for(int q=0;q<orbit_num_p;q++)
					{
						y[p][q]=0;
					}
				}
				for(int veci=0;veci<dim_p;veci++)
				{
					int ii=Q_coll_stru_tab_p[Qp_index].i[veci];
					int jj=Q_coll_stru_tab_p[Qp_index].j[veci];
					y[ii][jj]+=vt[veci+k*dim_p];
				}
				cout<<"Qp="<<endl;
				for(int p=0;p<orbit_num_p;p++)
				{
					for(int q=0;q<orbit_num_p;q++)
					{
						cout<<y[p][q]<<' ';
						QQ_qp_np[QQ_num_np][p][q]=y[p][q];
					}
					cout<<endl;
				}
				cout<<endl;
				for(int p=0;p<orbit_num_p;p++)
				{
					for(int q=p;q<orbit_num_p;q++)
					{
						double delta=y[p][q];
						if((orbit_p[p].j-orbit_p[q].j)%4==0)
						{
							delta-=y[q][p];
						}
						else
						{
							delta+=y[q][p];
						}
						
						if(fabs(delta)>1e-7)
						{
							isher=false;
						}
						delta=y[p][q];
						if((orbit_p[p].j+orbit_n[q].j)%4==0)
						{
							delta-=y[q][p];
						}
						else
						{
							delta+=y[q][p];
						}
						
						if(fabs(delta)>1e-7)
						{
							isantiher=false;
						}
					}
				}
				if(isher)
				{
					QQ_her_np[QQ_num_np]=1;
				}
				if(isantiher)
				{
					QQ_her_np[QQ_num_np]=-1;
				}
				if((!isher)&&(!isantiher))
				{
					QQ_her_np[QQ_num_np]=0;
				}
				QQ_num_np++;
			}
			delete [] u;
			delete [] vt;
			delete [] work;
		}
		delete[] h;
		delete[] eig;
	}
	for(int i=0;i<tb_num;i++)
	{
		if(tb_index[i][0]<orbit_num_n+1&&tb_index[i][1]>orbit_num_n+1)
		{
			if(fabs(tb_value[i]-tb_value_test[i])>1e-3)
			{
				for(int k=0;k<5;k++)
				{
					cout<<"error: "<<tb_index[i][k]<<' ';
				}
				cout<<tb_value[i]<<' '<<tb_value_test[i]<<endl;
			}
		}
	}
}

void dig_P_eig_n()
{
	int target[5];
	size_t target_index;
	char a = 'N';
	char b = 'U';
	int info = 0;
	double *h;
	double *eig;
	double *work;
	int lwork;
	int dim;
	int dim2;
	int inc = 1;
	P_eig_num_n = 0;
	int dim_max = -1;
	for (int i = 0; i < P_coll_stru_num_n; i++)
	{
		if (P_coll_stru_tab_n[i].dim > dim_max)
		{
			dim_max = P_coll_stru_tab_n[i].dim;
		}
	}
	h = new double[dim_max * dim_max];
	eig = new double[dim_max];
	work = new double[3 * dim_max];
	for (int i = 0; i < P_coll_stru_num_n; i++)
	{
		dim = P_coll_stru_tab_n[i].dim;
		dim2 = dim * dim;
		memset(h, 0, sizeof(double) * dim2);
		lwork = dim * 3;
		target[4] = P_coll_stru_tab_n[i].J / 2;
		for (int k = 0; k < dim; k++)
		{
			for (int l = k; l < dim; l++)
			{
				target[0] = P_coll_stru_tab_n[i].i[k] + 1;
				target[1] = P_coll_stru_tab_n[i].j[k] + 1;
				target[2] = P_coll_stru_tab_n[i].i[l] + 1;
				target[3] = P_coll_stru_tab_n[i].j[l] + 1;

				if (array_search(tb_num_n, (void **)tb_index_n, (void *)target, &target_index, compar))
				{
					h[k + l * dim] += tb_P_value_n[target_index];
				}
			}
		}
		bool flag_dia_able = false;
		for (int k = 0; k < dim2; k++)
		{
			if (fabs(h[k]) > 0)
			{
				flag_dia_able = true;
				break;
			}
		}
		if (flag_dia_able)
		{
			if (flag_print)
			{
				cout << "J=" << P_coll_stru_tab_n[i].J << endl;
				for (int k = 0; k < dim; k++)
				{
					for (int l = 0; l < dim; l++)
					{
						cout << h[k + l * dim] << ' ';
					}
					cout << endl;
				}
			}
			dsyev_(&a, &b, &dim, h, &dim, eig, work, &lwork, &info);
			if (flag_print)
			{
				cout << "eig=";
				for (int k = 0; k < dim; k++)
				{
					cout << eig[k] << ' ';
				}
				cout << endl;
			}
			memcpy(P_eig_n + P_eig_num_n, eig, sizeof(double) * dim);
			P_eig_num_n += dim;
		}
	}
	delete[] h;
	delete[] eig;
	delete[] work;
}


double *M_Q2P_n;
void Q2P_init_n()
{
	M_Q2P_n = new double[tb_num_n * Fk_num_n];
	int target[5];
	size_t target_index_for_Fk;
	size_t target_index_for_tb;
	int sign;
	double delta;
	double *Fk_I = new double[Fk_num_n];
	memset(Fk_I, 0, sizeof(double) * Fk_num_n);
	for (int Fk_index = 0; Fk_index < Fk_num_n; Fk_index++)
	{
		Fk_I[Fk_index] = 1;
		if (Fk_index > 0)
		{
			Fk_I[Fk_index - 1] = 0;
		}
		memset(M_Q2P_n + Fk_index * tb_num_n, 0, sizeof(double) * tb_num_n);
		for (int i = 0; i < orbit_num_n; i++)
		{
			for (int j = 0; j < orbit_num_n; j++)
			{
				for (int kappa = abs(orbit_n[i].j - orbit_n[j].j); kappa <= orbit_n[i].j + orbit_n[j].j; kappa += 2)
				{
					for (int k = 0; k < orbit_num_n; k++)
					{
						for (int l = 0; l < orbit_num_n; l++)
						{
							if (kappa < abs(orbit_n[k].j - orbit_n[l].j) || kappa > orbit_n[k].j + orbit_n[l].j)
							{
								continue;
							}
							if ((orbit_n[i].l + orbit_n[j].l + orbit_n[k].l + orbit_n[l].l) % 2 != 0)
							{
								continue;
							}
							target[0] = i + 1;
							target[1] = j + 1;
							target[2] = k + 1;
							target[3] = l + 1;
							target[4] = kappa / 2;
							sign = 1;
							if (target[0] > target[1])
							{
								if (abs(orbit_n[i].j - orbit_n[j].j) % 4 != 0)
								{
									sign *= -1;
								}

								swap(target, target + 1);
							}
							if (target[2] > target[3])
							{
								if (abs(orbit_n[k].j - orbit_n[l].j) % 4 != 0)
								{
									sign *= -1;
								}
								swap(target + 2, target + 3);
							}
							if (target[0] > target[2] || (target[0] == target[2] && target[1] > target[3]))
							{
								swap(target, target + 2);
								swap(target + 1, target + 3);
							}
							if (array_search(Fk_num_n, (void **)Fk_index_n, (void *)target, &target_index_for_Fk, compar))
							{
								if(fabs(Fk_I[target_index_for_Fk])<1e-15)
								{
									continue;
								}
								target[0] = i + 1;
								target[1] = k + 1;
								target[2] = j + 1;
								target[3] = l + 1;
								if (target[0] > target[1])
								{
									swap(target, target + 1);
								}
								if (target[2] > target[3])
								{
									swap(target + 2, target + 3);
								}
								if (target[0] > target[2] || (target[0] == target[2] && target[1] > target[3]))
								{
									continue;
								}
								for (int J = abs(orbit_n[i].j - orbit_n[k].j); J <= orbit_n[i].j + orbit_n[k].j; J += 2)
								{
									delta = sign * Fk_I[target_index_for_Fk];
									if (J < abs(orbit_n[j].j - orbit_n[l].j) || J > orbit_n[j].j + orbit_n[l].j)
									{
										continue;
									}
									delta *= (kappa + 1) * sixj_out(orbit_n[i].j, orbit_n[j].j, kappa, orbit_n[l].j, orbit_n[k].j, J);
									if ((orbit_n[j].j + orbit_n[k].j + J) % 4 != 0)
									{
										delta *= -1;
									}
									target[0] = i + 1;
									target[1] = k + 1;
									target[2] = j + 1;
									target[3] = l + 1;
									target[4] = J / 2;
									if (target[0] > target[1])
									{
										if (abs(orbit_n[i].j + orbit_n[k].j + J) % 4 == 0)
										{
											delta *= -1;
										}
										swap(target, target + 1);
									}
									if (target[2] > target[3])
									{
										if (abs(orbit_n[j].j + orbit_n[l].j + J) % 4 == 0)
										{
											delta *= -1;
										}
										swap(target + 2, target + 3);
									}
									if (i == k)
									{
										delta *= sqrt(2);
									}
									if (j == l)
									{
										delta *= sqrt(2);
									}
									if (array_search(tb_num_n, (void **)tb_index_n, (void *)target, &target_index_for_tb, compar))
									{
										(M_Q2P_n+ Fk_index * tb_num_n)[target_index_for_tb] += delta;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	delete[] Fk_I;
}

void Q2P_n()
{
	char t='N';
	int M=tb_num_n;
	int N=Fk_num_n;
	double alpha=1;
	int inc=1;
	double beta=0;
	dgemv_(&t,&M,&N,&alpha,M_Q2P_n,&M,Fk_value_n,&inc,&beta,tb_Q_value_n,&inc);
}

struct PQ_eig_vec *Q_eig_vec_n;
double *Q_eig_n;
int Q_eig_num_n;
int Q_eig_vec_num_n;
int *Q_eig_order_n;

int compar_Q_eig_vec_n(const void *p1, const void *p2)
{
	double eig1 = fabs(Q_eig_vec_n[*((int *)p1)].eig);
	double eig2 = fabs(Q_eig_vec_n[*((int *)p2)].eig);
	if (eig1 > eig2)
	{
		return -1;
	}
	else if (eig1 < eig2)
	{
		return 1;
	}
	return 0;
}

int compar_eig(const void *p1, const void *p2)
{
	double eig1 = fabs(*((double *)p1));
	double eig2 = fabs(*((double *)p2));
	if (eig1 > eig2)
	{
		return -1;
	}
	else if (eig1 < eig2)
	{
		return 1;
	}
	return 0;
}

void dig_Q_init_n()
{
	int dim;
	Q_eig_vec_num_n = 0;
	for (int i = 0; i < Q_coll_stru_num_n; i++)
	{
		dim = Q_coll_stru_tab_n[i].dim;
		Q_eig_vec_num_n += dim;
	}
	Q_eig_vec_n = new struct PQ_eig_vec[Q_eig_vec_num_n];
	Q_eig_order_n = new int[Q_eig_vec_num_n];
	Q_eig_n = new double[Q_eig_vec_num_n];
	Q_eig_num_n = Q_eig_vec_num_n;
	Q_eig_vec_num_n = 0;
	for (int i = 0; i < Q_coll_stru_num_n; i++)
	{
		dim = Q_coll_stru_tab_n[i].dim;
		for (int k = 0; k < dim; k++)
		{
			Q_eig_vec_n[Q_eig_vec_num_n].vec = new double[dim];
			Q_eig_vec_n[Q_eig_vec_num_n].coll_index = i;
			Q_eig_order_n[Q_eig_vec_num_n] = Q_eig_vec_num_n;
			Q_eig_vec_num_n++;
		}
	}
}
void dig_Q_eig_vec_n()
{
	int target[5];
	size_t target_index;
	char a = 'V';
	char b = 'U';
	int lwork;
	double *work;
	double *h;
	double *eig;
	int info = 0;
	int dim;
	int dim2;
	int inc = 1;
	int ii, jj, kk, ll, kappa, pi, J;
	int sign;
	int temp_i;
	double delta;
	Q_eig_vec_num_n = 0;
	for (int i = 0; i < Q_coll_stru_num_n; i++)
	{
		dim = Q_coll_stru_tab_n[i].dim;
		dim2 = dim * dim;
		h = new double[dim2];
		eig = new double[dim];
		lwork = 3 * dim;
		work = new double[lwork];
		memset(h, 0, sizeof(double) * dim2);
		kappa = Q_coll_stru_tab_n[i].J;
		for (int k = 0; k < dim; k++)
		{
			for (int l = k; l < dim; l++)
			{
				ii = Q_coll_stru_tab_n[i].i[k];
				jj = Q_coll_stru_tab_n[i].j[k];
				kk = Q_coll_stru_tab_n[i].i[l];
				ll = Q_coll_stru_tab_n[i].j[l];
				sign = 1;
				if (ii > jj)
				{
					swap(&ii, &jj);
					if ((orbit_n[ii].j - orbit_n[jj].j) % 4 != 0)
					{
						sign *= -1;
					}
				}
				if (kk > ll)
				{
					swap(&kk, &ll);
					if ((orbit_n[kk].j - orbit_n[ll].j) % 4 != 0)
					{
						sign *= -1;
					}
				}
				if (ii > kk || (ii == kk && jj > ll))
				{
					swap(&ii, &kk);
					swap(&jj, &ll);
				}
				target[0] = ii + 1;
				target[1] = jj + 1;
				target[2] = kk + 1;
				target[3] = ll + 1;
				target[4] = kappa / 2;
				if (array_search(Fk_num_n, (void **)Fk_index_n, (void *)target, &target_index, compar))
				{
					h[k + l * dim] += sign * Fk_value_n[target_index];
				}
				//cout<<h[k+l*dim]<<' ';
			}
			//cout<<endl;
		}
		bool flag_dia_able = false;
		for (int k = 0; k < dim2; k++)
		{
			if (fabs(h[k]) > 0)
			{
				flag_dia_able = true;
				break;
			}
		}
		if (flag_dia_able)
		{
			if (flag_print)
			{
				cout << "kappa=" << Q_coll_stru_tab_n[i].J << endl;
				for (int k = 0; k < dim; k++)
				{
					for (int l = 0; l < dim; l++)
					{
						cout << h[k + l * dim] << ' ';
					}
					cout << endl;
				}
			}
			dsyev_(&a, &b, &dim, h, &dim, eig, work, &lwork, &info);
			if (flag_print)
			{
				cout << "eig=";
				for (int k = 0; k < dim; k++)
				{
					cout << eig[k] << ' ';
				}
				cout << endl;
			}
			for (int k = 0; k < dim; k++)
			{
				Q_eig_vec_n[Q_eig_vec_num_n].eig = eig[k];
				Q_eig_vec_n[Q_eig_vec_num_n].coll_index = i;
				memcpy(Q_eig_vec_n[Q_eig_vec_num_n].vec, h + k * dim, sizeof(double) * dim);
				Q_eig_vec_num_n++;
			}
		}
		else
		{
			for (int k = 0; k < dim; k++)
			{
				Q_eig_vec_n[Q_eig_vec_num_n].eig = 0;
				Q_eig_vec_n[P_eig_vec_num_n].coll_index = i;
				memset(Q_eig_vec_n[Q_eig_vec_num_n].vec, 0, sizeof(double) * dim);
				Q_eig_vec_num_n++;
			}
		}

		delete[] work;
		delete[] h;
		delete[] eig;
	}
}
void dig_Q_eig_n()
{
	int target[5];
	size_t target_index;
	char a = 'N';
	char b = 'U';
	int lwork;
	double *work;
	double *h;
	double *eig;
	int info = 0;
	int dim;
	int dim2;
	int inc = 1;
	int ii, jj, kk, ll, kappa, pi, J;
	int sign;
	int temp_i;
	double delta;
	Q_eig_vec_num_n = 0;
	int dim_max = -1;
	for (int i = 0; i < Q_coll_stru_num_n; i++)
	{
		if (Q_coll_stru_tab_n[i].dim > dim_max)
		{
			dim_max = Q_coll_stru_tab_n[i].dim;
		}
	}
	h = new double[dim_max * dim_max];
	eig = new double[dim_max];
	work = new double[3 * dim_max];
	for (int i = 0; i < Q_coll_stru_num_n; i++)
	{
		dim = Q_coll_stru_tab_n[i].dim;
		dim2 = dim * dim;
		lwork = 3 * dim;
		memset(h, 0, sizeof(double) * dim2);
		kappa = Q_coll_stru_tab_n[i].J;
		for (int k = 0; k < dim; k++)
		{
			for (int l = k; l < dim; l++)
			{
				ii = Q_coll_stru_tab_n[i].i[k];
				jj = Q_coll_stru_tab_n[i].j[k];
				kk = Q_coll_stru_tab_n[i].i[l];
				ll = Q_coll_stru_tab_n[i].j[l];
				sign = 1;
				if (ii > jj)
				{
					swap(&ii, &jj);
					if ((orbit_n[ii].j - orbit_n[jj].j) % 4 != 0)
					{
						sign *= -1;
					}
				}
				if (kk > ll)
				{
					swap(&kk, &ll);
					if ((orbit_n[kk].j - orbit_n[ll].j) % 4 != 0)
					{
						sign *= -1;
					}
				}
				if (ii > kk || (ii == kk && jj > ll))
				{
					swap(&ii, &kk);
					swap(&jj, &ll);
				}
				target[0] = ii + 1;
				target[1] = jj + 1;
				target[2] = kk + 1;
				target[3] = ll + 1;
				target[4] = kappa / 2;
				if (array_search(Fk_num_n, (void **)Fk_index_n, (void *)target, &target_index, compar))
				{
					h[k + l * dim] += sign * Fk_value_n[target_index];
				}
				//cout<<h[k+l*dim]<<' ';
			}
			//cout<<endl;
		}
		bool flag_dia_able = false;
		for (int k = 0; k < dim2; k++)
		{
			if (fabs(h[k]) > 0)
			{
				flag_dia_able = true;
				break;
			}
		}
		if (flag_dia_able)
		{
			if (flag_print)
			{
				cout << "kappa=" << Q_coll_stru_tab_n[i].J << endl;
				for (int k = 0; k < dim; k++)
				{
					for (int l = 0; l < dim; l++)
					{
						cout << h[k + l * dim] << ' ';
					}
					cout << endl;
				}
			}
			dsyev_(&a, &b, &dim, h, &dim, eig, work, &lwork, &info);
			if (flag_print)
			{
				cout << "eig=";
				for (int k = 0; k < dim; k++)
				{
					cout << eig[k] << ' ';
				}
				cout << endl;
			}
			memcpy(Q_eig_n + Q_eig_vec_num_n, eig, sizeof(double) * dim);
			Q_eig_vec_num_n += dim;
		}
	}
	delete[] work;
	delete[] h;
	delete[] eig;
}
int P_cut;
int Q_cut;
void P2tb_n(int P_cut)
{
	memset(tb_P_value_n, 0, sizeof(double) * tb_num_n);
	int dim;
	int coll_index;
	int ii, jj, kk, ll, J;
	int target[5];
	size_t target_index;
	int index_to_vec;
	for (int i = 0; i < P_cut; i++)
	{
		index_to_vec = P_eig_order_n[i];
		coll_index = P_eig_vec_n[index_to_vec].coll_index;
		dim = P_coll_stru_tab_n[coll_index].dim;
		target[4] = P_coll_stru_tab_n[coll_index].J / 2;
		for (int p = 0; p < dim; p++)
		{
			for (int q = p; q < dim; q++)
			{
				target[0] = P_coll_stru_tab_n[coll_index].i[p] + 1;
				target[1] = P_coll_stru_tab_n[coll_index].j[p] + 1;
				target[2] = P_coll_stru_tab_n[coll_index].i[q] + 1;
				target[3] = P_coll_stru_tab_n[coll_index].j[q] + 1;
				if (array_search(tb_num_n, (void **)tb_index_n, (void *)target, &target_index, compar))
				{
					tb_P_value_n[target_index] += P_eig_vec_n[index_to_vec].eig * P_eig_vec_n[index_to_vec].vec[p] * P_eig_vec_n[index_to_vec].vec[q];
				}
			}
		}
	}
}
void Q2tb_n(int Q_cut)
{
	memset(Fk_value_n, 0, sizeof(double) * Fk_num_n);
	int dim;
	int coll_index;
	int ii, jj, kk, ll, J;
	int target[5];
	size_t target_index;
	int index_to_vec;
	for (int i = 0; i < Q_cut; i++)
	{
		index_to_vec = Q_eig_order_n[i];
		coll_index = Q_eig_vec_n[index_to_vec].coll_index;
		dim = Q_coll_stru_tab_n[coll_index].dim;
		target[4] = Q_coll_stru_tab_n[coll_index].J / 2;
		for (int p = 0; p < dim; p++)
		{
			if (Q_coll_stru_tab_n[coll_index].i[p] > Q_coll_stru_tab_n[coll_index].j[p])
			{
				continue;
			}
			for (int q = p; q < dim; q++)
			{
				if (Q_coll_stru_tab_n[coll_index].i[q] > Q_coll_stru_tab_n[coll_index].j[q])
				{
					continue;
				}
				target[0] = Q_coll_stru_tab_n[coll_index].i[p] + 1;
				target[1] = Q_coll_stru_tab_n[coll_index].j[p] + 1;
				target[2] = Q_coll_stru_tab_n[coll_index].i[q] + 1;
				target[3] = Q_coll_stru_tab_n[coll_index].j[q] + 1;
				if (array_search(Fk_num_n, (void **)Fk_index_n, (void *)target, &target_index, compar))
				{
					Fk_value_n[target_index] += Q_eig_vec_n[index_to_vec].eig * Q_eig_vec_n[index_to_vec].vec[p] * Q_eig_vec_n[index_to_vec].vec[q];
				}
			}
		}
	}
	Q2P_n();
}

double *Qm_sp_n;
double *Qm_sp_p;
double ***PPm_n;
double ***QQm_n;
double ***PPm_p;
double ***QQm_p;
double ***QQm_qn_np;
double ***QQm_qp_np;

void npa_ham_to_PQm()
{
	Qm_sp_n=new double [C_dim_n*C_dim_n];
	sp_to_Qm(Qm_sp_n,ob_n,C_dim_n,C_n);
	PPm_n=new double **[PP_num_n];
	for(int i=0;i<PP_num_n;i++)
	{
		PPm_n[i]=new double *[PP_J_n[i]+1];
		for(int k=0;k<PP_J_n[i]+1;k++)
		{
			PPm_n[i][k]=new double [C_dim_n*C_dim_n];
			P_old_npa_to_Pm(PPm_n[i][k],PP_y_n[i],PP_J_n[i],-PP_J_n[i]+k*2,C_dim_n,P_dim_n,P_stru_n,orbit_num_n,orbit_n);
		}
	}
	QQm_n=new double **[QQ_num_n];
	for(int i=0;i<QQ_num_n;i++)
	{
		QQm_n[i]=new double *[QQ_kappa_n[i]+1];
		for(int k=0;k<QQ_kappa_n[i]+1;k++)
		{
			QQm_n[i][k]=new double [C_dim_n*C_dim_n];
			Q_old_npa_to_Qm(QQm_n[i][k],QQ_q_n[i],QQ_kappa_n[i],-QQ_kappa_n[i]+k*2,C_dim_n,Q_dim_n,Q_stru_n,orbit_num_n,orbit_n);
		}
	}
	Qm_sp_p=new double [C_dim_p*C_dim_p];
	sp_to_Qm(Qm_sp_p,ob_p,C_dim_p,C_p);
	PPm_p=new double **[PP_num_p];
	for(int i=0;i<PP_num_p;i++)
	{
		PPm_p[i]=new double *[PP_J_p[i]+1];
		for(int k=0;k<PP_J_p[i]+1;k++)
		{
			PPm_p[i][k]=new double [C_dim_p*C_dim_p];
			P_old_npa_to_Pm(PPm_p[i][k],PP_y_p[i],PP_J_p[i],-PP_J_p[i]+k*2,C_dim_p,P_dim_p,P_stru_p,orbit_num_p,orbit_p);
		}
	}
	QQm_p=new double **[QQ_num_p];
	for(int i=0;i<QQ_num_p;i++)
	{
		QQm_p[i]=new double *[QQ_kappa_p[i]+1];
		for(int k=0;k<QQ_kappa_p[i]+1;k++)
		{
			QQm_p[i][k]=new double [C_dim_p*C_dim_p];
			Q_old_npa_to_Qm(QQm_p[i][k],QQ_q_p[i],QQ_kappa_p[i],-QQ_kappa_p[i]+k*2,C_dim_p,Q_dim_p,Q_stru_p,orbit_num_p,orbit_p);
		}
	}
	QQm_qn_np=new double **[QQ_num_np];
	QQm_qp_np=new double **[QQ_num_np];
	for(int i=0;i<QQ_num_np;i++)
	{
		QQm_qn_np[i]=new double *[QQ_kappa_np[i]+1];
		QQm_qp_np[i]=new double *[QQ_kappa_np[i]+1];
		for(int k=0;k<QQ_kappa_np[i]+1;k++)
		{
			QQm_qn_np[i][k]=new double [C_dim_n*C_dim_n];
			Q_old_npa_to_Qm(QQm_qn_np[i][k],QQ_qn_np[i],QQ_kappa_np[i],-QQ_kappa_np[i]+k*2,C_dim_n,Q_dim_n,Q_stru_n,orbit_num_n,orbit_n);
			QQm_qp_np[i][k]=new double [C_dim_p*C_dim_p];
			Q_old_npa_to_Qm(QQm_qp_np[i][k],QQ_qp_np[i],QQ_kappa_np[i],-QQ_kappa_np[i]+k*2,C_dim_p,Q_dim_p,Q_stru_p,orbit_num_p,orbit_p);
		}
	}
}

/*
void PQ_eig_vec_out_n(const gsl_vector *v, void *par)
{
	for (int i = 0; i < Fk_num_n; i++)
	{
		Fk_value_n[i] = gsl_vector_get(v, i);
	}
	dig_Q_eig_vec_n();
	Q2P_n();
	double alpha = -1.0;
	int inc = 1;
	memcpy(tb_P_value_n, tb_value_n, sizeof(double) * tb_num_n);
	daxpy(&tb_num_n, &alpha, tb_Q_value_n, &inc, tb_P_value_n, &inc);
	dig_P_eig_vec_n();
	qsort(P_eig_order_n, P_eig_vec_num_n, sizeof(int), compar_P_eig_vec_n);
	qsort(Q_eig_order_n, Q_eig_vec_num_n, sizeof(int), compar_Q_eig_vec_n);
	double *temp = new double[tb_num_n];
	double rms;
	int index;
	int coll_index;
	Q2tb_n(Q_cut);
	P2tb_n(P_cut);
	memcpy(temp, tb_Q_value_n, sizeof(double) * tb_num_n);
	alpha = 1;
	daxpy_(&tb_num_n, &alpha, tb_P_value_n, &inc, temp, &inc);
	alpha = -1;
	daxpy_(&tb_num_n, &alpha, tb_value_n, &inc, temp, &inc);
	rms = dnrm2_(&tb_num_n, temp, &inc);
	rms *= rms;
	rms /= tb_num_n;
	rms = (sqrt(rms));
	cout<<"the final Fk= ";
	for(int i=0;i<Fk_num_n;i++)
	{
		cout<<Fk_value_n[i]<<',';
	}
	cout<<endl;
	cout<<"the resultant rms= "<<rms<<" with P_cut= "<<P_cut<<" Q_cut= "<<Q_cut<<endl;
	cout<<"so you get the final H with below structur:"<<endl;
	cout<<"P="<<endl;
	double **yq=new double * [orbit_num_n];
	for(int i=0;i<orbit_num_n;i++)
	{
		yq[i]=new double [orbit_num_n];
	}
	for(int i=0;i<P_cut;i++)
	{
		for(int k=0;k<orbit_num_n;k++)
		{
			memset(yq[k],0,sizeof(double)*orbit_num_n);
		}
		index=P_eig_order_n[i];
		coll_index=P_eig_vec_n[index].coll_index;
		cout<<i<<"th : "<<P_coll_stru_tab_n[coll_index].J<<' '<<P_coll_stru_tab_n[coll_index].pi<<' '<<P_eig_vec_n[index].eig<<endl;
		int *Ai=P_coll_stru_tab_n[coll_index].i;
		int *Aj=P_coll_stru_tab_n[coll_index].j;
		double *vec=P_eig_vec_n[index].vec;
		double eig=P_eig_vec_n[index].eig;
		int dim=P_coll_stru_tab_n[coll_index].dim;
		int J=P_coll_stru_tab_n[coll_index].J;
		cout<<"vec= ";
		for(int k=0;k<P_coll_stru_tab_n[coll_index].dim;k++)
		{
			cout<<P_eig_vec_n[index].vec[k]<<' ';
		}
		cout<<endl;
		for(int k=0;k<dim;k++)
		{
			if(Ai[k]==Aj[k])
			{
				yq[Ai[k]][Aj[k]]+=vec[k]/sqrt(2);
			}
			else
			{
				yq[Ai[k]][Aj[k]]+=vec[k]/2;
				if((orbit_n[Ai[k]].j+orbit_n[Aj[k]].j+J)%4==0)
				{
					yq[Aj[k]][Ai[k]]+=vec[k]/2;
				}
				else
				{
					yq[Aj[k]][Ai[k]]-=vec[k]/2;
				}
			}
		}
		cout<<"y[][]="<<endl;
		for(int p=0;p<orbit_num_n;p++)
		{
			for(int q=0;q<orbit_num_n;q++)
			{
				cout<<yq[p][q]<<',';
			}
			cout<<endl;
		}
	}
	int Q0_index=JP_in_PQ_coll_stru_tab(0,0,Q_coll_stru_num_n,Q_coll_stru_tab_n);
	int Q0_dim=Q_coll_stru_tab_n[Q0_index].dim;
	int *Q0_i=Q_coll_stru_tab_n[Q0_index].i;
	int *Q0_j=Q_coll_stru_tab_n[Q0_index].j;
	double *Q0_vec=new double [Q0_dim];
	double delta;
	memset(Q0_vec,0,sizeof(double)*Q0_dim);
	cout<<"Q="<<endl;
	for(int i=0;i<Q_cut;i++)
	{
		for(int k=0;k<orbit_num_n;k++)
		{
			memset(yq[k],0,sizeof(double)*orbit_num_n);
		}
		index=Q_eig_order_n[i];
		coll_index=Q_eig_vec_n[index].coll_index;
		cout<<i<<"th : "<<Q_coll_stru_tab_n[coll_index].J<<' '<<P_coll_stru_tab_n[coll_index].pi<<' '<<Q_eig_vec_n[index].eig<<endl;
		int *Qi=Q_coll_stru_tab_n[coll_index].i;
		int *Qj=Q_coll_stru_tab_n[coll_index].j;
		int dim=Q_coll_stru_tab_n[coll_index].dim;
		double *Qvec=Q_eig_vec_n[index].vec;
		double eig=Q_eig_vec_n[index].eig;
		int kappa=Q_coll_stru_tab_n[coll_index].J;
		cout<<"vec= ";
		for(int k=0;k<Q_coll_stru_tab_n[coll_index].dim;k++)
		{
			cout<<Q_eig_vec_n[index].vec[k]<<' ';
		}
		cout<<endl;
		for(int k=0;k<dim;k++)
		{
			yq[Qi[k]][Qj[k]]+=Qvec[k];
		}
		cout<<"q[][]="<<endl;
		for(int p=0;p<orbit_num_n;p++)
		{
			for(int q=0;q<orbit_num_n;q++)
			{
				cout<<yq[p][q]<<',';
			}
			cout<<endl;
		}
		for(int p=0;p<Q_coll_stru_tab_n[coll_index].dim;p++)
		{
			for(int q=0;q<Q_coll_stru_tab_n[coll_index].dim;q++)
			{
				if(orbit_n[Qi[p]].j==orbit_n[Qj[q]].j&&Qj[p]==Qi[q])
				{
					Q0_index=ij_in_Q_coll_stru(Qi[p],Qj[q],0,0,orbit_num_n,orbit_n);
					delta=eig*Qvec[p]*Qvec[q]*(kappa+1)/sqrt(orbit_n[Qi[p]].j+1);
					if((orbit_n[Qi[p]].j+orbit_n[Qj[p]].j)%4!=0)
					{
						delta=-delta;
					}
					Q0_vec[Q0_index]-=delta;
				}
			}
		}
	}
	for(int k=0;k<orbit_num_n;k++)
	{
		memset(yq[k],0,sizeof(double)*orbit_num_n);
	}
	for(int k=0;k<Q0_dim;k++)
	{
		yq[Q0_i[k]][Q0_j[k]]+=Q0_vec[k];
	}
	cout<<"Q0+="<<endl;
	for(int p=0;p<orbit_num_n;p++)
	{
		for(int q=0;q<orbit_num_n;q++)
		{
			if((orbit_n[p].l+orbit_n[q].l)%2==0)
			{
				cout<<yq[p][q]<<',';
			}
			else
			{
				cout<<0<<',';
			}
			
		}
		cout<<endl;
	}
	cout<<"Q0-="<<endl;
	for(int p=0;p<orbit_num_n;p++)
	{
		for(int q=0;q<orbit_num_n;q++)
		{
			if((orbit_n[p].l+orbit_n[q].l)%2!=0)
			{
				cout<<yq[p][q]<<',';
			}
			else
			{
				cout<<0<<',';
			}
			
		}
		cout<<endl;
	}
	exit(0);
	delete[] temp;
	delete [] Q0_vec;
}


double PQ_eig_vec_fit_n(const gsl_vector *v, void *par)
{
	for (int i = 0; i < Fk_num_n; i++)
	{
		Fk_value_n[i] = gsl_vector_get(v, i);
	}
	dig_Q_eig_vec_n();
	Q2P_n();
	double alpha = -1.0;
	int inc = 1;
	memcpy(tb_P_value_n, tb_value_n, sizeof(double) * tb_num_n);
	daxpy(&tb_num_n, &alpha, tb_Q_value_n, &inc, tb_P_value_n, &inc);
	dig_P_eig_vec_n();
	qsort(P_eig_order_n, P_eig_vec_num_n, sizeof(int), compar_P_eig_vec_n);
	qsort(Q_eig_order_n, Q_eig_vec_num_n, sizeof(int), compar_Q_eig_vec_n);
	double *temp = new double[tb_num_n];
	double rms;
	double rms_min=1000;
	
	int P_eig_order_index;
	int Q_eig_order_index;
	int index;
	int coll_index;
	//for(Q_cut=0;Q_cut<Q_eig_vec_num_n;Q_cut++)
	{
		Q2tb_n(Q_cut);
		//for(P_cut=0;P_cut<P_eig_vec_num_n;P_cut++)
		{
			P2tb_n(P_cut);
			memcpy(temp, tb_Q_value_n, sizeof(double) * tb_num_n);
			alpha = 1;
			daxpy_(&tb_num_n, &alpha, tb_P_value_n, &inc, temp, &inc);
			alpha = -1;
			daxpy_(&tb_num_n, &alpha, tb_value_n, &inc, temp, &inc);
			rms = dnrm2_(&tb_num_n, temp, &inc);
			rms *= rms;
			rms /= tb_num_n;
			rms = sqrt(rms);
			// if(rms<rms_min)
			// {
			// 	rms_min=rms;
			// 	//cout<<Q_cut<<' '<<P_cut<<' '<<pow(10,rms-Q_cut-P_cut)<<endl;
			// }
		}
	}
	
	delete[] temp;
	return rms;
}
gsl_vector *v_in_df;
int df_index;

double f_for_df_PQ_eig_vec_n(double x, void *par)
{
	double v_bak = gsl_vector_get(v_in_df, df_index);
	gsl_vector_set(v_in_df, df_index, x);
	double res = PQ_eig_vec_fit_n(v_in_df, NULL);
	gsl_vector_set(v_in_df, df_index, v_bak);
	return res;
}

void df_PQ_eig_vec_n(const gsl_vector *v, void *par, gsl_vector *df)
{
	gsl_function F;
	double result, abserr;

	F.function = &f_for_df_PQ_eig_vec_n;
	F.params = NULL;
	gsl_vector_memcpy(v_in_df, v);
	for (int i = 0; i < Fk_num_n; i++)
	{
		df_index = i;
		gsl_deriv_central(&F, gsl_vector_get(v, i), 1e-8, &result, &abserr);
		gsl_vector_set(df, i, result);
	}
}

void fdf_PQ_eig_vec_n(const gsl_vector *x, void *par,
						  double *f, gsl_vector *df)
{
	*f = PQ_eig_vec_fit_n(x, NULL);
	df_PQ_eig_vec_n(x, NULL, df);
}

gsl_vector *ss;
gsl_multimin_fminimizer *s_no_df;
const gsl_multimin_fminimizer_type *T_no_df;
double min_PQ_eig_vec_no_df_n(gsl_vector *x)
{
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status;
	double size;


	gsl_vector_set_all(ss, 0.1);

	minex_func.n = Fk_num_n;
	minex_func.f = PQ_eig_vec_fit_n;
	minex_func.params = NULL;

	gsl_multimin_fminimizer_set(s_no_df, &minex_func, x, ss);
	clock_t start, finish;
	start = clock();
	bool isfirst=true;
	double f_bak;
	do
	{
		iter++;
		status = gsl_multimin_fminimizer_iterate(s_no_df);

		if (status)
			break;

		size = gsl_multimin_fminimizer_size(s_no_df);
		status = gsl_multimin_test_size(size, 0.00001);
		if(fabs(f_bak-s_no_df->fval)<1e-5&&!isfirst)
		{
			cout<<"do not go any further"<<endl;
			break;
		}
		isfirst=false;
		f_bak=s_no_df->fval;
		if (status == GSL_SUCCESS)
		{
			printf("converged to minimum at\n");
		}
		finish = clock();
		if (double(finish - start) / CLOCKS_PER_SEC > 5)
		{
			start = clock();
			printf("%d %d %d f() = %f size = %f\n", P_cut,Q_cut,iter, s_no_df->fval, size);
			// if(f_rec-s->f>0.001)
			// {
			// 	f_rec=s->f;
			// }
			// else
			// {
			// 	cout<<"it's too slow to stick around"<<endl;
			// 	break;
			// }
		}

	} while (status == GSL_CONTINUE && iter < 100000000000);
	gsl_vector_memcpy(x, s_no_df->x);
	return s_no_df->fval;
}
gsl_multimin_fdfminimizer *s;
const gsl_multimin_fdfminimizer_type *T;
double min_PQ_eig_vec_n(gsl_vector *x)
{

	gsl_multimin_function_fdf my_func;

	my_func.n = Fk_num_n;
	my_func.f = PQ_eig_vec_fit_n;
	my_func.df = df_PQ_eig_vec_n;
	my_func.fdf = fdf_PQ_eig_vec_n;
	my_func.params = NULL;

	size_t iter = 0;
	int status;
	double size;

	clock_t start, finish;
	start = clock();
	double f_bak=1000;
	int same_num=0;


	//for(eig_cut=1;eig_cut<=Fk_num_n;eig_cut++)
	{
		//gsl_vector_set_zero(x);
		gsl_multimin_fdfminimizer_set(s, &my_func, x, 0.1, 1e-1);
		//f_bak=s->f;
		do
		{
			iter++;
			status = gsl_multimin_fdfminimizer_iterate(s);

			if (status)
				break;

			status = gsl_multimin_test_gradient(s->gradient, 1e-5);
			if(fabs(f_bak-s->f)<1e-8&&same_num>1)
			{
				cout<<"do not go any further"<<endl;
				break;
			}
			if(fabs(f_bak-s->f)<1e-8)
			{
				same_num++;
			}
			f_bak=s->f;
			if (status == GSL_SUCCESS)
			{
				printf("converged to minimum at\n");
				printf("%ld f() = %f \n", iter, s->f);
			}
			finish = clock();
			//if(double(finish-start)/CLOCKS_PER_SEC>10)
			{
				start = clock();
				printf("%d %d %ld f()= %f \n", P_cut,Q_cut,iter, s->f);
				// for(int i=0;i<Fk_num_n;i++)
				// {
				// 	cout<<gsl_vector_get(s->gradient,i)<<' ';
				// }
				// cout<<endl;
			}

		} while (status == GSL_CONTINUE && iter < 1000);
	}
	gsl_vector_memcpy(x, s->x);
	return s->f;
}
*/