int Q_tran_num[2];
int *Q_tran_her[2];
int *Q_M_tran_index[2];
int tran_num;
struct tran_struct
{
	int J_l;
	int pi_l;
	int i_l;
	int Q_norp;
	int Q_i;
	int J_r;
	int pi_r;
	int i_r;
};
vector<tran_struct> tran;
void tran_input_read()
{
	FILE *fp;
	char filename[80] = "tran_input.dat";
	double temp_d;
	if ((fp = fopen(filename, "r")) == NULL)
	{
		cout << "open file error" << endl;
	}
	else
	{
		fscanf(fp, "%d", &tran_num);
		for (int i = 0; i < tran_num; i++)
		{
			tran_struct temp;
			fscanf(fp, "%d %d %d %d %d %d %d %d", &(temp.J_l), &(temp.pi_l), &(temp.i_l), &(temp.Q_norp), &(temp.Q_i), &(temp.J_r), &(temp.pi_r), &(temp.i_r));
			temp.Q_norp-=1;
			temp.Q_i-=1;
			tran.emplace_back(temp);
		}
		for (int norp = 0; norp < 2; norp++)
		{
			fscanf(fp, "%d", &(Q_tran_num[norp]));
			Q_tran_her[norp] = new int[Q_tran_num[norp]];
			Q_M_tran_index[norp] = new int[Q_tran_num[norp]];
			for (int i = 0; i < Q_tran_num[norp]; i++)
			{
				PQ_Jpi_struct Q_Jpi_temp;
				fscanf(fp, "%d %d %d", &(Q_Jpi_temp.J), &(Q_Jpi_temp.pi), Q_tran_her[norp] + i);
				Q_Jpi_temp.yq = new double[sp_nlj_dim2[norp]];
				for (int p = 0; p < sp_nlj_dim2[norp]; p++)
				{
					fscanf(fp, "%lf", Q_Jpi_temp.yq + p);
				}
				Q_Jpi[norp].emplace_back(Q_Jpi_temp);
				Q_M_tran_index[norp][i] = QJpi_to_QM(Q_Jpi[norp].size() - 1, 0, norp) + INT_MAX / 2;
			}
		}
		fclose(fp);
	}
}