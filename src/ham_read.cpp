int PP_num[2];
int QQ_num[2];
int QQ_np_num;

vector<double> ob[2];
vector<int> Q_ob_index[2];
struct PP_QQ_ham_struct
{
	int J;
	int pi;
	double par;
	vector<int> PQ_index;
};
vector<PP_QQ_ham_struct> PP_ham[2];

vector<PP_QQ_ham_struct> QQ_ham[2];

struct QQ_np_ham_struct
{
	int J;
	int pi;
	double par;
	int her;
	vector<int > Q_index[2];
};
vector<QQ_np_ham_struct> QQ_np_ham;
void ham_read()
{
	FILE *fp;
	char filename[80]="ham.dat";
	double temp_d;
	if((fp=fopen(filename,"r"))==NULL)
	{
		cout<<"open file error"<<endl;
	}
	else
	{
		for(int norp=0;norp<2;norp++)
		{
			for(int i=0;i<sp_nlj[norp].size();i++)
			{
				fscanf(fp, "%lf", &temp_d);
				ob[norp].emplace_back(temp_d);
			}
			for(int k=0;k<sp_nlj[norp].size();k++)
			{
				double *q=new double [sp_nljm_dim2[norp]];
				memset(q,0,sizeof(double)*sp_nljm_dim2[norp]);

				for(int i=0;i<sp_nljm[norp].size();i++)
				{
					for(int j=0;j<sp_nljm[norp].size();j++)
					{
						if(sp_nljm[norp][i].m==sp_nljm[norp][j].m&&sp_nljm[norp][i].nlj==sp_nljm[norp][j].nlj&&sp_nljm[norp][i].nlj==k)
						{
							q[i+j*sp_nljm[norp].size()]=1;
						}
					}
				}
				PQ_M_struct temp;
				temp.J=0;
				temp.M=0;
				temp.pi=0;
				temp.PQ_Jpi_index=-1;
				temp.yq=q;
				temp.yq_mM=q;
				Q_M[norp].emplace_back(temp);
				Q_ob_index[norp].emplace_back(Q_M[norp].size()-1+INT_MAX/2);
				// Q_ob_index[norp][k]=Q_M[norp].size()-1+INT_MAX/2;
			}
		}		
		for(int norp=0;norp<2;norp++)
		{
			fscanf(fp,"%d",&(PP_num[norp]));
			for(int i=0;i<PP_num[norp];i++)
			{
				PP_QQ_ham_struct PP_ham_temp;
				PQ_Jpi_struct P_Jpi_temp;
				fscanf(fp,"%lf %d %d",&(PP_ham_temp.par),&(PP_ham_temp.J),&(PP_ham_temp.pi));
				P_Jpi_temp.J=PP_ham_temp.J;
				P_Jpi_temp.pi=PP_ham_temp.pi;
				P_Jpi_temp.yq=new double [sp_nlj_dim2[norp]];
				for(int p=0;p<sp_nlj_dim2[norp];p++)
				{
					fscanf(fp,"%lf", P_Jpi_temp.yq+p);
				}
				P_Jpi[norp].emplace_back(P_Jpi_temp);
				for(int M=0;M<=P_Jpi_temp.J;M+=2)
				{
					int P_index=PJpi_to_PM(P_Jpi[norp].size()-1,M,norp);
					PP_ham_temp.PQ_index.emplace_back(P_index);
				}
				PP_ham[norp].emplace_back(PP_ham_temp);
			}
			fscanf(fp,"%d",&(QQ_num[norp]));
			for(int i=0;i<QQ_num[norp];i++)
			{
				PP_QQ_ham_struct QQ_ham_temp;
				PQ_Jpi_struct Q_Jpi_temp;
				fscanf(fp,"%lf %d %d",&(QQ_ham_temp.par),&(QQ_ham_temp.J),&(QQ_ham_temp.pi));
				Q_Jpi_temp.J=QQ_ham_temp.J;
				Q_Jpi_temp.pi=QQ_ham_temp.pi;
				Q_Jpi_temp.yq=new double [sp_nlj_dim2[norp]];
				for(int p=0;p<sp_nlj_dim2[norp];p++)
				{
					fscanf(fp,"%lf", Q_Jpi_temp.yq+p);
				}
				Q_Jpi[norp].emplace_back(Q_Jpi_temp);
				for(int M=0;M<=Q_Jpi_temp.J;M+=2)
				{
					int Q_index=QJpi_to_QM(Q_Jpi[norp].size()-1,M,norp);
					QQ_ham_temp.PQ_index.emplace_back(Q_index+INT_MAX/2);
				}
				QQ_ham[norp].emplace_back(QQ_ham_temp);
			}
		}
		
		fscanf(fp,"%d",&QQ_np_num);
		for(int i=0;i<QQ_np_num;i++)
		{
			QQ_np_ham_struct QQ_np_temp;
			fscanf(fp,"%lf %d %d %d\n",&(QQ_np_temp.par),&(QQ_np_temp.J),&(QQ_np_temp.pi),&(QQ_np_temp.her));
			PQ_Jpi_struct Q1_Jpi_temp,Q2_Jpi_temp;
			Q1_Jpi_temp.J=QQ_np_temp.J;
			Q2_Jpi_temp.J=QQ_np_temp.J;
			Q1_Jpi_temp.pi=QQ_np_temp.pi;
			Q2_Jpi_temp.pi=QQ_np_temp.pi;
			Q1_Jpi_temp.yq=new double [sp_nlj_dim2[0]];
			Q2_Jpi_temp.yq=new double [sp_nlj_dim2[1]];
			for(int p=0;p<sp_nlj_dim2[0];p++)
			{
				fscanf(fp,"%lf",Q1_Jpi_temp.yq+p);
			}
			for(int p=0;p<sp_nlj_dim2[1];p++)
			{
				fscanf(fp,"%lf",Q2_Jpi_temp.yq+p);
			}
			Q_Jpi[0].emplace_back(Q1_Jpi_temp);
			Q_Jpi[1].emplace_back(Q2_Jpi_temp);
			for(int M=0;M<=QQ_np_temp.J;M+=2)
			{
				int Q1_index=QJpi_to_QM(Q_Jpi[0].size()-1,M,0);
				QQ_np_temp.Q_index[0].emplace_back(Q1_index+INT_MAX/2);
				int Q2_index=QJpi_to_QM(Q_Jpi[1].size()-1,M,1);
				QQ_np_temp.Q_index[1].emplace_back(Q2_index+INT_MAX/2);
			}
			QQ_np_ham.emplace_back(QQ_np_temp);
		}
		fclose(fp);
	}
}