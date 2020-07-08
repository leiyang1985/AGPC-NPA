struct block_mat_struct
{
	int pi;
	int M;
	int dim;
	double *ove;
	double *vec;
	double *JpJm;
	vector<double *> spe;
	vector<double *> PP;
	vector<double *> QQ;
	vector<double *> PP_in_basis;
	vector<double *> C_in_basis;
	double *JpJm_unnor;
	vector<double *> spe_unnor;
	vector<double *> PP_unnor;
	vector<double *> QQ_unnor;
	vector<double *> PP_in_basis_unnor;
	vector<double *> C_in_basis_unnor;
};
vector<block_mat_struct> block_mat[2][2];
int block_store_num[2][2] = {{0, 0}, {0, 0}};
vector<block_mat_struct> block_store[2];
int *block_store_index_in_mat[2];
bool *block_mat_ove_evalued[2][2];
bool *block_mat_ham_evalued[2][2];
void block_mat_init()
{
	for (int norp = 0; norp < 2; norp++)
	{
		for (int pi = 0; pi < 2; pi++)
		{
			block_mat[norp][pi].resize(basis_block[norp][pi].size());
			block_mat_ove_evalued[norp][pi] = new bool[basis_block[norp][pi].size()];
			block_mat_ham_evalued[norp][pi] = new bool[basis_block[norp][pi].size()];
			for (int i = 0; i < basis_block[norp][pi].size(); i++)
			{
				block_mat_ove_evalued[norp][pi][i] = false;
				block_mat_ham_evalued[norp][pi][i] = false;
			}
		}
	}
}
void block_mat_read()
{
	char filename[] = "block_mat_0_0.bin";
	for (int norp = 0; norp < 2; norp++)
	{
		for (int pi = 0; pi < 2; pi++)
		{
			filename[10] = '0' + norp;
			filename[12] = '0' + pi;
			FILE *fp;
			if((fp= fopen(filename, "rb"))==NULL)
			{
				block_store_num[norp][pi]=0;
				continue;
			}
			fread(block_store_num[norp] + pi, sizeof(int), 1, fp);
			for (int i = 0; i < block_store_num[norp][pi]; i++)
			{
				
				int M, dim;
				fread(&(M), sizeof(int), 1, fp);
				fread(&(dim), sizeof(int), 1, fp);
				int count_int=2;
				int count_double=0;
				bool Mpi_need = false;
				int Mpi_index;
				for (int k = 0; k < basis_block[norp][pi].size(); k++)
				{
					if (basis_block[norp][pi][k].M == M)
					{
						Mpi_need = true;
						Mpi_index = k;
						block_mat[norp][pi][k].M = M;
						block_mat[norp][pi][k].dim = dim;
						block_mat_ove_evalued[norp][pi][k] = true;
						block_mat_ham_evalued[norp][pi][k] = true;
						break;
					}
				}
				int dim2=dim*(dim+1)/2;
				if (Mpi_need)
				{
					if (dim < basis_block[norp][pi][Mpi_index].r.size())
					{
						basis_block[norp][pi][Mpi_index].r.resize(dim);
						for (int k = 0; k < dim; k++)
						{
							fread(&(basis_block[norp][pi][Mpi_index].r[k][0]), sizeof(int),N[norp] / 2, fp);
							count_int+=N[norp]/2;
						}
						if(M==0)
						{
							fread(Rev_p_M0[norp][pi],sizeof(int),dim,fp);
							count_int+=dim;
						}
					}
					else
					{
						fseek(fp, dim * (N[norp] / 2) * sizeof(int), SEEK_CUR);
						count_int+=dim * (N[norp] / 2);
						if(M==0)
						{
							fseek(fp, dim * sizeof(int), SEEK_CUR);
							count_int+=dim;
						}
					}
					
					
					block_mat[norp][pi][Mpi_index].ove = new double[dim2];
					fread(block_mat[norp][pi][Mpi_index].ove, sizeof(double), dim2, fp);
					block_mat[norp][pi][Mpi_index].vec = new double[dim2];
					fread(block_mat[norp][pi][Mpi_index].vec, sizeof(double), dim2, fp);
					block_mat[norp][pi][Mpi_index].JpJm = new double[dim2];
					fread(block_mat[norp][pi][Mpi_index].JpJm, sizeof(double), dim2, fp);
					count_double+=dim2*3;
					for (int k = 0; k < sp_nlj[norp].size(); k++)
					{
						double *temp_dp = new double[dim2];
						fread(temp_dp, sizeof(double), dim2, fp);
						block_mat[norp][pi][Mpi_index].spe.emplace_back(temp_dp);
						count_double+=dim2;
					}
					for (int k = 0; k < PP_num[norp]; k++)
					{
						double *temp_dp = new double[dim2];
						fread(temp_dp, sizeof(double), dim2, fp);
						block_mat[norp][pi][Mpi_index].PP.emplace_back(temp_dp);
						count_double+=dim2;
					}
					for (int k = 0; k < QQ_num[norp]; k++)
					{
						double *temp_dp = new double[dim2];
						fread(temp_dp, sizeof(double), dim2, fp);
						block_mat[norp][pi][Mpi_index].QQ.emplace_back(temp_dp);
						count_double+=dim2;
					}
					for(int k=0;k<P_Jpi_in_basis_num[norp];k++)
					{
						double *temp_dp=new double [dim2];
						fread(temp_dp,sizeof(double),dim2,fp);
						block_mat[norp][pi][Mpi_index].PP_in_basis.emplace_back(temp_dp);
						count_double+=dim2;
					}
					block_mat[norp][pi][Mpi_index].JpJm_unnor = new double[dim2];
					fread(block_mat[norp][pi][Mpi_index].JpJm_unnor, sizeof(double), dim2, fp);
					count_double+=dim2;
					for (int k = 0; k < sp_nlj[norp].size(); k++)
					{
						double *temp_dp = new double[dim2];
						fread(temp_dp, sizeof(double), dim2, fp);
						block_mat[norp][pi][Mpi_index].spe_unnor.emplace_back(temp_dp);
						count_double+=dim2;
					}
					for (int k = 0; k < PP_num[norp]; k++)
					{
						double *temp_dp = new double[dim2];
						fread(temp_dp, sizeof(double), dim2, fp);
						block_mat[norp][pi][Mpi_index].PP_unnor.emplace_back(temp_dp);
						count_double+=dim2;
					}
					for (int k = 0; k < QQ_num[norp]; k++)
					{
						double *temp_dp = new double[dim2];
						fread(temp_dp, sizeof(double), dim2, fp);
						block_mat[norp][pi][Mpi_index].QQ_unnor.emplace_back(temp_dp);
						count_double+=dim2;
					}
					for(int k=0;k<P_Jpi_in_basis_num[norp];k++)
					{
						double *temp_dp=new double [dim2];
						fread(temp_dp,sizeof(double),dim2,fp);
						block_mat[norp][pi][Mpi_index].PP_in_basis_unnor.emplace_back(temp_dp);
						count_double+=dim2;
					}
					if(N[norp]%2==1)
					{
						for(int k=0;k<sp_nlj[norp].size();k++)
						{
							double *temp_dp=new double [dim2];
							fread(temp_dp,sizeof(double),dim2,fp);
							block_mat[norp][pi][Mpi_index].C_in_basis.emplace_back(temp_dp);
							count_double+=dim2;
						}
						for(int k=0;k<sp_nlj[norp].size();k++)
						{
							double *temp_dp=new double [dim2];
							fread(temp_dp,sizeof(double),dim2,fp);
							block_mat[norp][pi][Mpi_index].C_in_basis_unnor.emplace_back(temp_dp);
							count_double+=dim2;
						}
					}
				}
				else
				{
					fseek(fp, dim * (N[norp] / 2) * sizeof(int), SEEK_CUR);
					count_int+=dim*(N[norp]/2);
					fseek(fp, (4 + 2*sp_nlj[norp].size() + 2*PP_num[norp] + 2*QQ_num[norp]+2*P_Jpi_in_basis_num[norp]) * dim2 * sizeof(double), SEEK_CUR);
					count_double+=(4 + 2*sp_nlj[norp].size() + 2*PP_num[norp] + 2*QQ_num[norp]+2*P_Jpi_in_basis_num[norp]) * dim2;
					if(N[norp]%2==1)
					{
						fseek(fp, 2*sp_nlj[norp].size() * dim2 * sizeof(double), SEEK_CUR);
						count_double+=2*sp_nlj[norp].size() * dim2;
					}
				}
				cout<<"reading block_mat: "<<norp<<' '<<pi<<' '<<i<<' '<<M<<' '<<dim<<' '<<count_int<<' '<<count_double<<endl;
			}
			fclose(fp);
		}
	}
}
void block_mat_write()
{
	char filename[] = "block_mat_0_0.bin";
	for (int norp = 0; norp < 2; norp++)
	{
		for (int pi = 0; pi < 2; pi++)
		{
			if (block_mat[norp][pi].size() <= block_store_num[norp][pi])
			{
				continue;
			}
			filename[10] = '0' + norp;
			filename[12] = '0' + pi;
			FILE *fp = fopen(filename, "wb");
			int temp_d = block_mat[norp][pi].size();
			fwrite(&temp_d, sizeof(int), 1, fp);
			for (int i = 0; i < temp_d; i++)
			{
				fwrite(&(basis_block[norp][pi][i].M), sizeof(int), 1, fp);
				int dim = basis_block[norp][pi][i].r.size();
				int count_int=2;
				int count_double=0;
				
				fwrite(&dim, sizeof(int), 1, fp);
				for (int k = 0; k < dim; k++)
				{
					fwrite(&(basis_block[norp][pi][i].r[k][0]), sizeof(int), N[norp] / 2, fp);
					count_int+=N[norp]/2;
				}
				if(basis_block[norp][pi][i].M==0)
				{
					fwrite(Rev_p_M0[norp][pi],sizeof(int),dim,fp);
					count_int+=dim;
				}
				int dim2 = dim * (dim + 1) / 2;
				fwrite(block_mat[norp][pi][i].ove, sizeof(double), dim2, fp);
				fwrite(block_mat[norp][pi][i].vec, sizeof(double), dim2, fp);
				fwrite(block_mat[norp][pi][i].JpJm, sizeof(double), dim2, fp);
				count_double+=dim2*3;
				for (int k = 0; k < sp_nlj[norp].size(); k++)
				{
					fwrite(block_mat[norp][pi][i].spe[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				for (int k = 0; k < PP_num[norp]; k++)
				{
					fwrite(block_mat[norp][pi][i].PP[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				for (int k = 0; k < QQ_num[norp]; k++)
				{
					fwrite(block_mat[norp][pi][i].QQ[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
				{
					fwrite(block_mat[norp][pi][i].PP_in_basis[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				fwrite(block_mat[norp][pi][i].JpJm_unnor, sizeof(double), dim2, fp);
				count_double+=dim2;
				for (int k = 0; k < sp_nlj[norp].size(); k++)
				{
					fwrite(block_mat[norp][pi][i].spe_unnor[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				for (int k = 0; k < PP_num[norp]; k++)
				{
					fwrite(block_mat[norp][pi][i].PP_unnor[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				for (int k = 0; k < QQ_num[norp]; k++)
				{
					fwrite(block_mat[norp][pi][i].QQ_unnor[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
				{
					fwrite(block_mat[norp][pi][i].PP_in_basis_unnor[k], sizeof(double), dim2, fp);
					count_double+=dim2;
				}
				if(N[norp]%2==1)
				{
					for (int k = 0; k < sp_nlj[norp].size(); k++)
					{
						fwrite(block_mat[norp][pi][i].C_in_basis[k], sizeof(double), dim2, fp);
						count_double+=dim2;
					}
					for (int k = 0; k < sp_nlj[norp].size(); k++)
					{
						fwrite(block_mat[norp][pi][i].C_in_basis_unnor[k], sizeof(double), dim2, fp);
						count_double+=dim2;
					}
				}
				cout<<"writing block_mat: "<<norp<<' '<<pi<<' '<<i<<' '<<basis_block[norp][pi][i].M<<' '<<dim<<' '<<count_int<<' '<<count_double<<endl;
			}
			fclose(fp);
		}
	}
}

void schmit_vec_out(double *ove, double *vec, int dim)
{
	for (int p = 0; p < dim; p++)
	{
		memset(vec + pack_k(0, p), 0, sizeof(double) * (p + 1));
		vec[pack_k(p, p)] = 1;
		for (int k = 0; k < p; k++)
		{
			double delta = 0;
			for (int i = 0; i < k + 1; i++)
			{
				delta += vec[pack_k(i, k)] * ove[pack_k(i, p)];
			}
			delta = -delta;
			int NN = k + 1;
			int inc = 1;
			daxpy_(&NN, &delta, vec + pack_k(0, k), &inc, vec + pack_k(0, p), &inc);
		}
		double sum2 = 0;
		for (int i = 0; i < p + 1; i++)
		{
			for (int j = i; j < p + 1; j++)
			{
				if (i != j)
				{
					sum2 += 2 * ove[pack_k(i, j)] * vec[pack_k(i, p)] * vec[pack_k(j, p)];
				}
				else
				{
					sum2 += ove[pack_k(i, j)] * vec[pack_k(i, p)] * vec[pack_k(j, p)];
				}
			}
		}
		if(sum2>0)
		{
			sum2 = 1 / sqrt(sum2);
		}
		else
		{
			sum2=0;
		}
		int NN = p + 1;
		int inc = 1;
		dscal_(&NN, &sum2, vec + pack_k(0, p), &inc);
	}
}
void basis_delete(int i, int &p, int norp, int pi)
{
	int M = basis_block[norp][pi][i].M;
	int dim=basis_block[norp][pi][i].r.size();
	int Rev_i = basis_block[norp][pi].size() - 1 - i;
	int Rev_p;
	if (M != 0)
	{
		Rev_p = p;
		basis_block[norp][pi][i].r.erase(basis_block[norp][pi][i].r.begin() + p);
		block_mat[norp][pi][i].dim--;
		basis_block[norp][pi][Rev_i].r.erase(basis_block[norp][pi][Rev_i].r.begin() + Rev_p);
		block_mat[norp][pi][Rev_i].dim--;
		p--;
	}
	else
	{
		Rev_p = Rev_p_M0[norp][pi][p];
		if (p > Rev_p)
		{
			basis_block[norp][pi][i].r.erase(basis_block[norp][pi][i].r.begin() + p);
			memmove(Rev_p_M0[norp][pi]+p,Rev_p_M0[norp][pi]+p+1,sizeof(int)*(dim-p-1));
			dim--;
			basis_block[norp][pi][Rev_i].r.erase(basis_block[norp][pi][Rev_i].r.begin() + Rev_p);
			memmove(Rev_p_M0[norp][pi]+Rev_p,Rev_p_M0[norp][pi]+Rev_p+1,sizeof(int)*(dim-Rev_p-1));
			dim--;
			for(int k=0;k<dim;k++)
			{
				if(Rev_p_M0[norp][pi][k]>p)
				{
					Rev_p_M0[norp][pi][k]-=2;
				}
				else if(Rev_p_M0[norp][pi][k]>Rev_p&&Rev_p_M0[norp][pi][k]<p)
				{
					Rev_p_M0[norp][pi][k]-=1;
				}
			}
			block_mat[norp][pi][i].dim = dim;
		}
		else if(p<Rev_p)
		{
			basis_block[norp][pi][Rev_i].r.erase(basis_block[norp][pi][Rev_i].r.begin() + Rev_p);
			memmove(Rev_p_M0[norp][pi]+Rev_p,Rev_p_M0[norp][pi]+Rev_p+1,sizeof(int)*(dim-Rev_p-1));
			dim--;
			basis_block[norp][pi][i].r.erase(basis_block[norp][pi][i].r.begin() + p);
			memmove(Rev_p_M0[norp][pi]+p,Rev_p_M0[norp][pi]+p+1,sizeof(int)*(dim-p-1));
			dim--;
			for(int k=0;k<dim;k++)
			{
				if(Rev_p_M0[norp][pi][k]>Rev_p)
				{
					Rev_p_M0[norp][pi][k]-=2;
				}
				else if(Rev_p_M0[norp][pi][k]<Rev_p&&Rev_p_M0[norp][pi][k]>p)
				{
					Rev_p_M0[norp][pi][k]-=1;
				}
			}
			block_mat[norp][pi][i].dim = dim;
		}
		else
		{
			basis_block[norp][pi][i].r.erase(basis_block[norp][pi][i].r.begin() + p);
			memmove(Rev_p_M0[norp][pi]+p,Rev_p_M0[norp][pi]+p+1,sizeof(int)*(dim-p-1));
			dim--;
			for(int k=0;k<dim;k++)
			{
				if(Rev_p_M0[norp][pi][k]>p)
				{
					Rev_p_M0[norp][pi][k]-=1;
				}
			}
			block_mat[norp][pi][i].dim=dim;
		}
		p--;
	}
}
void basis_filter(int norp, int pi)
{
	PQ_mat_struct sr[(N[norp] / 2 +N[norp]%2)* 2];
	for (int i = 0; i < basis_block[norp][pi].size(); i++)
	{
		if (block_mat_ove_evalued[norp][pi][i])
		{
			continue;
		}
		int M, dim, dim2;
		M = basis_block[norp][pi][i].M;
		dim = basis_block[norp][pi][i].r.size();
		block_mat[norp][pi][i].M = M;
		block_mat[norp][pi][i].dim = dim;
		dim2 = dim * (dim + 1) / 2;
		block_mat[norp][pi][i].ove = new double[dim2];
		block_mat[norp][pi][i].vec = new double[dim2];
		int Rev_i = basis_block[norp][pi].size() - 1 - i;
		int Rev_p;
		int Rev_q;
		double delta;
		if (block_mat_ove_evalued[norp][pi][Rev_i])
		{
			for (int p = 0; p < basis_block[norp][pi][i].r.size(); p++)
			{
				if (basis_block[norp][pi][i].M != 0)
				{
					Rev_p = p;
				}
				else
				{
					Rev_p = Rev_p_M0[norp][pi][p];
				}
				for (int q = p; q < basis_block[norp][pi][i].r.size(); q++)
				{
					if (basis_block[norp][pi][i].M != 0)
					{
						Rev_q = q;
					}
					else
					{
						Rev_q = Rev_p_M0[norp][pi][q];
					}
					bool swaped = false;
					if (Rev_p > Rev_q)
					{
						swaped = true;
						swap(&Rev_p, &Rev_q);
					}
					int phase = 0;
					for (int k = 0; k < N[norp] / 2; k++)
					{
						phase += P_M[norp][abs(basis_block[norp][pi][i].r[p][k])].J;
						phase += P_M[norp][abs(basis_block[norp][pi][i].r[q][k])].J;
					}
					if(N[norp]%2!=0)
					{
						phase+=sp_nljm[norp][abs(basis_block[norp][pi][i].r[p][N[norp]/2])-1].j;
						phase+=sp_nljm[norp][abs(basis_block[norp][pi][i].r[q][N[norp]/2])-1].j;
					}
					phase-=2*M;
					if (phase % 4 == 0)
					{
						phase = 1;
					}
					else
					{
						phase = -1;
					}
					block_mat[norp][pi][i].ove[pack_k(p, q)] = phase * block_mat[norp][pi][Rev_i].ove[pack_k(Rev_p, Rev_q)];
					if (swaped)
					{
						swap(&Rev_p, &Rev_q);
					}
				}
			}
			schmit_vec_out(block_mat[norp][pi][i].ove, block_mat[norp][pi][i].vec, dim);
			block_mat_ove_evalued[norp][pi][i] = true;
		}
		else
		{
			cout << norp << ' ' << pi << " valence space's " << i << "th block has M=" << basis_block[norp][pi][i].M << " with dim=" << basis_block[norp][pi][i].r.size() << " in basis_filter"<< endl;
			for (int p = 0; p < basis_block[norp][pi][i].r.size(); p++)
			{
				for (int k = 0; k < N[norp] / 2; k++)
				{
					// if (basis_block[norp][pi][i].r[p][k] > 0 || P_M[norp][pi][basis_block[norp][pi][i].r[p][k]].M)
					if (basis_block[norp][pi][i].r[p][k] > 0)
					{
						sr[k] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[p][k]].J, P_M[norp][basis_block[norp][pi][i].r[p][k]].pi, P_M[norp][basis_block[norp][pi][i].r[p][k]].M, P_M[norp][basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
					}
					else
					{
						sr[k] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[p][k]].J, P_M[norp][-basis_block[norp][pi][i].r[p][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[p][k]].M, P_M[norp][-basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
					}
					if (basis_block[norp][pi][i].r[p][k] > 0)
					{
						sr[k + N[norp] / 2+N[norp]%2] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[p][k]].J, P_M[norp][basis_block[norp][pi][i].r[p][k]].pi, P_M[norp][basis_block[norp][pi][i].r[p][k]].M, P_M[norp][basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
					}
					else
					{
						sr[k + N[norp] / 2+N[norp]%2] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[p][k]].J, P_M[norp][-basis_block[norp][pi][i].r[p][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[p][k]].M, P_M[norp][-basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
					}
				}
				if(N[norp]%2!=0)
				{
					if(basis_block[norp][pi][i].r[p][N[norp]/2]>0)
					{
						int sp_index=basis_block[norp][pi][i].r[p][N[norp]/2]-1;
						sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,0);
					}
					else
					{
						int sp_index=-basis_block[norp][pi][i].r[p][N[norp]/2]-1;
						sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),0);
					}
					if(basis_block[norp][pi][i].r[p][N[norp]/2]>0)
					{
						int sp_index=basis_block[norp][pi][i].r[p][N[norp]/2]-1;
						sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,1);
					}
					else
					{
						int sp_index=-basis_block[norp][pi][i].r[p][N[norp]/2]-1;
						sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),1);
					}
				}
				if(N[norp]%2==0)
				{
					delta = ove_cal(sr, N[norp] / 2, norp);
				}
				else
				{
					delta = ove_o_cal(sr, norp);
				}
				cout << "ove=" << delta << endl;
				if (delta > 1e-14)
				{
					block_mat[norp][pi][i].ove[pack_k(p, p)] = delta;
					for (int q = 0; q < p; q++)
					{
						for (int k = 0; k < N[norp] / 2; k++)
						{
							if (basis_block[norp][pi][i].r[p][k] > 0)
							{
								sr[k] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[p][k]].J, P_M[norp][basis_block[norp][pi][i].r[p][k]].pi, P_M[norp][basis_block[norp][pi][i].r[p][k]].M, P_M[norp][basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
							}
							else
							{
								sr[k] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[p][k]].J, P_M[norp][-basis_block[norp][pi][i].r[p][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[p][k]].M, P_M[norp][-basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
							}
							if (basis_block[norp][pi][i].r[q][k] > 0)
							{
								sr[k + N[norp] / 2+N[norp]%2] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[q][k]].J, P_M[norp][basis_block[norp][pi][i].r[q][k]].pi, P_M[norp][basis_block[norp][pi][i].r[q][k]].M, P_M[norp][basis_block[norp][pi][i].r[q][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[q][k]);
							}
							else
							{
								sr[k + N[norp] / 2+N[norp]%2] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[q][k]].J, P_M[norp][-basis_block[norp][pi][i].r[q][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[q][k]].M, P_M[norp][-basis_block[norp][pi][i].r[q][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[q][k]);
							}
						}
						if(N[norp]%2==0)
						{
							sr_sort(sr,N[norp]/2);
							block_mat[norp][pi][i].ove[pack_k(q, p)] = ove_cal(sr, N[norp] / 2, norp);
						}
						else
						{
							if(basis_block[norp][pi][i].r[p][N[norp]/2]>0)
							{
								int sp_index=basis_block[norp][pi][i].r[p][N[norp]/2]-1;
								sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,0);
							}
							else
							{
								int sp_index=-basis_block[norp][pi][i].r[p][N[norp]/2]-1;
								sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),0);
							}
							if(basis_block[norp][pi][i].r[q][N[norp]/2]>0)
							{
								int sp_index=basis_block[norp][pi][i].r[q][N[norp]/2]-1;
								sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,1);
							}
							else
							{
								int sp_index=-basis_block[norp][pi][i].r[q][N[norp]/2]-1;
								sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),1);
							}
							sr_o_sort(sr,N[norp]/2);
							block_mat[norp][pi][i].ove[pack_k(q, p)] = ove_o_cal(sr, norp);
						}
						cout << block_mat[norp][pi][i].ove[pack_k(q, p)] << ' ';
					}
					cout << endl;
					memset(block_mat[norp][pi][i].vec + pack_k(0, p), 0, sizeof(double) * (p));
					block_mat[norp][pi][i].vec[pack_k(p, p)] = 1;
					for (int k = 0; k < p; k++)
					{
						double delta = 0;
						for (int ii = 0; ii < k + 1; ii++)
						{
							delta += block_mat[norp][pi][i].vec[pack_k(ii, k)] * block_mat[norp][pi][i].ove[pack_k(ii, p)];
						}
						delta = -delta;
						int NN = k + 1;
						int inc = 1;
						daxpy_(&NN, &delta, block_mat[norp][pi][i].vec + pack_k(0, k), &inc, block_mat[norp][pi][i].vec + pack_k(0, p), &inc);
					}
					double sum2 = 0;
					for (int ii = 0; ii < p + 1; ii++)
					{
						for (int jj = ii; jj < p + 1; jj++)
						{
							if(ii==jj)
							{
								sum2 += block_mat[norp][pi][i].ove[pack_k(ii, jj)] * block_mat[norp][pi][i].vec[pack_k(ii, p)] * block_mat[norp][pi][i].vec[pack_k(jj, p)];
							}
							if (ii != jj)
							{
								sum2 += 2*block_mat[norp][pi][i].ove[pack_k(ii, jj)] * block_mat[norp][pi][i].vec[pack_k(jj, p)] * block_mat[norp][pi][i].vec[pack_k(ii, p)];
							}
						}
					}
					if (sum2 > 1e-9)
					{
						sum2 = 1 / sqrt(sum2);
						int NN = p + 1;
						int inc = 1;
						dscal_(&NN, &sum2, block_mat[norp][pi][i].vec + pack_k(0, p), &inc);
					}
					else
					{
						cout << "delete this basis" << endl;
						basis_delete(i, p, norp, pi);
					}
				}
				else
				{
					cout << "delete this basis" << endl;
					basis_delete(i, p, norp, pi);
				}
			}
			block_mat_ove_evalued[norp][pi][i] = true;
		}
	}
}
void ham_mat_ort_nor(double *mat_out,double *mat, double *vec, int dim)
{
	int dim2 = dim * (dim + 1) / 2;
	memset(mat_out, 0, sizeof(double) * dim2);
	for (int p = 0; p < dim; p++)
	{
		for (int q = p; q < dim; q++)
		{
			for (int ii = 0; ii < dim; ii++)
			{
				for (int jj = ii; jj < dim; jj++)
				{
					if (ii <= p && jj <= q)
					{
						mat_out[pack_k(p, q)] += mat[pack_k(ii, jj)] * vec[pack_k(ii, p)] * vec[pack_k(jj, q)];
					}
					if (ii != jj && jj <= p && ii <= q)
					{
						mat_out[pack_k(p, q)] += mat[pack_k(ii, jj)] * vec[pack_k(jj, p)] * vec[pack_k(ii, q)];
					}
				}
			}
		}
	}
}

void block_mat_out(int norp, int pi)
{
	PQ_mat_struct sr[(N[norp] / 2 +N[norp]%2)* 2];
	for (int i = 0; i < basis_block[norp][pi].size(); i++)
	{
		if (block_mat_ham_evalued[norp][pi][i])
		{
			continue;
		}
		int M, dim, dim2;
		M = basis_block[norp][pi][i].M;
		dim = basis_block[norp][pi][i].r.size();
		dim2 = dim * (dim + 1) / 2;
		block_mat[norp][pi][i].JpJm = new double[dim2];
		memset(block_mat[norp][pi][i].JpJm,0,sizeof(double)*dim2);
		block_mat[norp][pi][i].spe.resize(sp_nlj[norp].size());
		for (int k = 0; k < sp_nlj[norp].size(); k++)
		{
			block_mat[norp][pi][i].spe[k] = new double[dim2];
			memset(block_mat[norp][pi][i].spe[k],0,sizeof(double)*dim2);
		}
		block_mat[norp][pi][i].PP.resize(PP_num[norp]);
		for (int k = 0; k < PP_num[norp]; k++)
		{
			block_mat[norp][pi][i].PP[k] = new double[dim2];
			memset(block_mat[norp][pi][i].PP[k],0,sizeof(double)*dim2);
		}
		block_mat[norp][pi][i].QQ.resize(PP_num[norp]);
		for (int k = 0; k < QQ_num[norp]; k++)
		{
			block_mat[norp][pi][i].QQ[k] = new double[dim2];
			memset(block_mat[norp][pi][i].QQ[k],0,sizeof(double)*dim2);
		}
		block_mat[norp][pi][i].PP_in_basis.resize(P_Jpi_in_basis_num[norp]);
		for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
		{
			block_mat[norp][pi][i].PP_in_basis[k] = new double[dim2];
			memset(block_mat[norp][pi][i].PP_in_basis[k],0,sizeof(double)*dim2);
		}
		block_mat[norp][pi][i].JpJm_unnor = new double[dim2];
		memset(block_mat[norp][pi][i].JpJm_unnor,0,sizeof(double)*dim2);
		block_mat[norp][pi][i].spe_unnor.resize(sp_nlj[norp].size());
		for (int k = 0; k < sp_nlj[norp].size(); k++)
		{
			block_mat[norp][pi][i].spe_unnor[k] = new double[dim2];
			memset(block_mat[norp][pi][i].spe_unnor[k],0,sizeof(double)*dim2);
		}
		block_mat[norp][pi][i].PP_unnor.resize(PP_num[norp]);
		for (int k = 0; k < PP_num[norp]; k++)
		{
			block_mat[norp][pi][i].PP_unnor[k] = new double[dim2];
			memset(block_mat[norp][pi][i].PP_unnor[k],0,sizeof(double)*dim2);
		}
		block_mat[norp][pi][i].QQ_unnor.resize(PP_num[norp]);
		for (int k = 0; k < QQ_num[norp]; k++)
		{
			block_mat[norp][pi][i].QQ_unnor[k] = new double[dim2];
			memset(block_mat[norp][pi][i].QQ_unnor[k],0,sizeof(double)*dim2);
		}
		block_mat[norp][pi][i].PP_in_basis_unnor.resize(P_Jpi_in_basis_num[norp]);
		for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
		{
			block_mat[norp][pi][i].PP_in_basis_unnor[k] = new double[dim2];
			memset(block_mat[norp][pi][i].PP_in_basis_unnor[k],0,sizeof(double)*dim2);
		}
		if(N[norp]%2==1)
		{
			block_mat[norp][pi][i].C_in_basis.resize(sp_nlj[norp].size());
			for (int k = 0; k < sp_nlj[norp].size(); k++)
			{
				block_mat[norp][pi][i].C_in_basis[k] = new double[dim2];
				memset(block_mat[norp][pi][i].C_in_basis[k],0,sizeof(double)*dim2);
			}
			block_mat[norp][pi][i].C_in_basis_unnor.resize(sp_nlj[norp].size());
			for (int k = 0; k < sp_nlj[norp].size(); k++)
			{
				block_mat[norp][pi][i].C_in_basis_unnor[k] = new double[dim2];
				memset(block_mat[norp][pi][i].C_in_basis_unnor[k],0,sizeof(double)*dim2);
			}
		}
		int Rev_i = basis_block[norp][pi].size() - 1 - i;
		if (block_mat_ham_evalued[norp][pi][Rev_i])
		{
			int Rev_p;
			int Rev_q;
			for (int p = 0; p < basis_block[norp][pi][i].r.size(); p++)
			{
				if (basis_block[norp][pi][i].M != 0)
				{
					Rev_p = p;
				}
				else
				{
					Rev_p = Rev_p_M0[norp][pi][p];
				}
				for (int q = p; q < basis_block[norp][pi][i].r.size(); q++)
				{
					if (basis_block[norp][pi][i].M != 0)
					{
						Rev_q = q;
					}
					else
					{
						Rev_q = Rev_p_M0[norp][pi][q];
					}
					bool swaped = false;
					if (Rev_p > Rev_q)
					{
						swap(&Rev_p, &Rev_q);
					}
					int phase = 0;
					for (int k = 0; k < N[norp] / 2; k++)
					{
						phase += P_M[norp][abs(basis_block[norp][pi][i].r[p][k])].J;
						phase += P_M[norp][abs(basis_block[norp][pi][i].r[q][k])].J;
					}
					if(N[norp]%2!=0)
					{
						phase+=sp_nljm[norp][abs(basis_block[norp][pi][i].r[p][N[norp]/2])-1].j;
						phase+=sp_nljm[norp][abs(basis_block[norp][pi][i].r[q][N[norp]/2])-1].j;
					}
					phase-=2*M;
					if (phase % 4 == 0)
					{
						phase = 1;
					}
					else
					{
						phase = -1;
					}
					block_mat[norp][pi][i].JpJm_unnor[pack_k(p, q)] = phase * block_mat[norp][pi][Rev_i].JpJm_unnor[pack_k(Rev_p, Rev_q)]+M*block_mat[norp][pi][i].ove[pack_k(p,q)];
					for (int k = 0; k < sp_nlj[norp].size(); k++)
					{
						block_mat[norp][pi][i].spe_unnor[k][pack_k(p, q)] = phase * block_mat[norp][pi][Rev_i].spe_unnor[k][pack_k(Rev_p, Rev_q)];
					}

					for (int k = 0; k < PP_num[norp]; k++)
					{
						block_mat[norp][pi][i].PP_unnor[k][pack_k(p, q)] = phase * block_mat[norp][pi][Rev_i].PP_unnor[k][pack_k(Rev_p, Rev_q)];
					}
					for (int k = 0; k < QQ_num[norp]; k++)
					{
						block_mat[norp][pi][i].QQ_unnor[k][pack_k(p, q)] = phase * block_mat[norp][pi][Rev_i].QQ_unnor[k][pack_k(Rev_p, Rev_q)];
					}
					for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
					{
						block_mat[norp][pi][i].PP_in_basis_unnor[k][pack_k(p, q)] = phase * block_mat[norp][pi][Rev_i].PP_in_basis_unnor[k][pack_k(Rev_p, Rev_q)];
					}
					if(N[norp]%2==1)
					{
						for (int k = 0; k < sp_nlj[norp].size(); k++)
						{
							block_mat[norp][pi][i].C_in_basis_unnor[k][pack_k(p, q)] = phase * block_mat[norp][pi][Rev_i].C_in_basis_unnor[k][pack_k(Rev_p, Rev_q)];
						}
					}
					if (swaped)
					{
						swap(&Rev_p, &Rev_q);
					}
				}
			}
		}
		else
		{
			cout << norp << ' ' << pi << " valence space's " << i << "th block has M=" << basis_block[norp][pi][i].M << " with dim=" << basis_block[norp][pi][i].r.size() << " in ham_mat_out"<<endl;
			cout<<"basis="<<endl;
			if(N[norp]%2==0)
			{
				for(int k=0;k<basis_block[norp][pi][i].r.size();k++)
				{
					for(int l=0;l<N[norp]/2;l++)
					{
						cout<<basis_block[norp][pi][i].r[k][l]<<' ';
					}
					cout<<endl;
					
				}
			}
			else
			{
				for(int k=0;k<basis_block[norp][pi][i].r.size();k++)
				{
					for(int l=0;l<N[norp]/2+1;l++)
					{
						cout<<basis_block[norp][pi][i].r[k][l]<<' ';
					}
					cout<<endl;
					
				}
			}
			
			cout<<"mat="<<endl;
			for(int p=0;p<basis_block[norp][pi][i].r.size();p++)
			{
				for(int q=p;q<basis_block[norp][pi][i].r.size();q++)
				{

					
					if(N[norp]%2==0)
					{
						for (int k = 0; k < N[norp] / 2; k++)
						{
							if (basis_block[norp][pi][i].r[p][k] > 0)
							{
								sr[k] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[p][k]].J, P_M[norp][basis_block[norp][pi][i].r[p][k]].pi, P_M[norp][basis_block[norp][pi][i].r[p][k]].M, P_M[norp][basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
							}
							else
							{
								sr[k] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[p][k]].J, P_M[norp][-basis_block[norp][pi][i].r[p][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[p][k]].M, P_M[norp][-basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
							}

							if (basis_block[norp][pi][i].r[q][k] > 0)
							{
								sr[k + N[norp] / 2] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[q][k]].J, P_M[norp][basis_block[norp][pi][i].r[q][k]].pi, P_M[norp][basis_block[norp][pi][i].r[q][k]].M, P_M[norp][basis_block[norp][pi][i].r[q][k]].PQ_Jpi_index,basis_block[norp][pi][i].r[q][k]);
							}
							else
							{
								sr[k + N[norp] / 2] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[q][k]].J, P_M[norp][-basis_block[norp][pi][i].r[q][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[q][k]].M, P_M[norp][-basis_block[norp][pi][i].r[q][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[q][k]);
							}
						}
						sr_sort(sr, N[norp] / 2);

						auto Q_Jp = PQ_mat_out(2, 0, 2, -1, J_pm_index[norp]);
						auto Q_Jm = PQ_mat_out(2, 0, -2, -1, -J_pm_index[norp]);
						block_mat[norp][pi][i].JpJm_unnor[pack_k(p, q)] = qq_cal(sr, Q_Jp, Q_Jm, N[norp] / 2, norp);
						
						PQ_mat_free(Q_Jp);
						PQ_mat_free(Q_Jm);

						for (int k = 0; k < sp_nlj[norp].size(); k++)
						{
							PQ_mat_struct Q0 = PQ_mat_out(0, 0, 0, -1, Q_ob_index[norp][k]);
							block_mat[norp][pi][i].spe_unnor[k][pack_k(p, q)] = q_cal(sr, Q0, N[norp] / 2, norp);
							PQ_mat_free(Q0);
						}
						for (int k = 0; k < PP_num[norp]; k++)
						{
							int J = PP_ham[norp][k].J;
							int pi_PP = PP_ham[norp][k].pi;
							int M_num = PP_ham[norp][k].PQ_index.size();
							for (int l = 0; l < M_num; l++)
							{
								int M = P_M[norp][PP_ham[norp][k].PQ_index[l]].M;
								int PQ_Jpi_index = P_M[norp][PP_ham[norp][k].PQ_index[l]].PQ_Jpi_index;
								auto P = PQ_mat_out(J, pi_PP, M, PQ_Jpi_index, PP_ham[norp][k].PQ_index[l]);
								double delta = pp_cal(sr, P, P, N[norp] / 2, norp);
								block_mat[norp][pi][i].PP_unnor[k][pack_k(p, q)] += delta;
								// block_mat_temp.PP[k][p+q*dim]+=pp_cal(sr,P,P,N[norp]/2,norp);
								PQ_mat_free(P);
								if (M != 0)
								{
									P = PQ_mat_out(J, pi_PP, -M, PQ_Jpi_index, -PP_ham[norp][k].PQ_index[l]);
									double delta = pp_cal(sr, P, P, N[norp] / 2, norp);

									block_mat[norp][pi][i].PP_unnor[k][pack_k(p, q)] += delta;
									// block_mat_temp.PP[k][p+q*dim]+=pp_cal(sr,P,P,N[norp]/2,norp);
									PQ_mat_free(P);
								}
							}
						}
						for (int k = 0; k < QQ_num[norp]; k++)
						{
							int J = QQ_ham[norp][k].J;
							int pi_QQ = QQ_ham[norp][k].pi;
							int M_num = QQ_ham[norp][k].PQ_index.size();
							for (int l = 0; l < M_num; l++)
							{
								int M = Q_M[norp][QQ_ham[norp][k].PQ_index[l] - INT_MAX / 2].M;
								int PQ_Jpi_index = Q_M[norp][QQ_ham[norp][k].PQ_index[l] - INT_MAX / 2].PQ_Jpi_index;
								auto Q = PQ_mat_out(J, pi_QQ, M, PQ_Jpi_index, QQ_ham[norp][k].PQ_index[l]);
								auto QT = PQ_mat_out(J, pi_QQ, -M, PQ_Jpi_index, -QQ_ham[norp][k].PQ_index[l]);
								double delta = qq_cal(sr, Q, QT, N[norp] / 2, norp);
								block_mat[norp][pi][i].QQ_unnor[k][pack_k(p, q)] += delta;
								// block_mat_temp.QQ[k][p+q*dim]+=qq_cal(sr,Q,QT,N[norp]/2,norp);
								if (M != 0)
								{
									delta = qq_cal(sr, QT, Q, N[norp] / 2, norp);
									block_mat[norp][pi][i].QQ_unnor[k][pack_k(p, q)] += delta;
									// block_mat_temp.QQ[k][p+q*dim]+=qq_cal(sr,QT,Q,N[norp]/2,norp);
								}
								PQ_mat_free(Q);
								PQ_mat_free(QT);
							}
						}
						for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
						{
							/*
							int J = P_Jpi[norp][k].J;
							int pi = P_Jpi[norp][k].pi;
							double pair_nor=0;
							for(int pp=0;pp<sp_nlj[norp].size();pp++)
							{
								for(int qq=pp;qq<sp_nlj[norp].size();qq++)
								{
									if(pp==qq)
									{
										pair_nor+=P_Jpi[norp][k].yq[pp+qq*sp_nlj[norp].size()]*P_Jpi[norp][k].yq[pp+qq*sp_nlj[norp].size()]*2;
									}
									else
									{
										pair_nor+=P_Jpi[norp][k].yq[pp+qq*sp_nlj[norp].size()]*P_Jpi[norp][k].yq[pp+qq*sp_nlj[norp].size()]*4;
									}
									
								}
							}
							for (int l = 0; l < P_in_basis[norp].size(); l++)
							{
								int P_M_index=P_in_basis[norp][l];
								if(P_M[norp][abs(P_M_index)].PQ_Jpi_index==k)
								{
									int M = P_M[norp][abs(P_M_index)].M;
									if(P_M_index<0)
									{
										M=-M;
									}
									auto P = PQ_mat_out(J, pi, M, k, P_M_index);
									double delta = pp_cal(sr, P, P, N[norp] / 2, norp)/pair_nor;
									block_mat[norp][pi][i].PP_in_basis_unnor[k][pack_k(p, q)] += delta;
									PQ_mat_free(P);
								}
								
							}
							*/
							int pair_k_num=0;
							for(int l=0;l<N[norp]/2;l++)
							{
								if(P_M[norp][abs(basis_block[norp][pi][i].r[p][l])].PQ_Jpi_index==k)
								{
									pair_k_num++;
								}
								if(P_M[norp][abs(basis_block[norp][pi][i].r[q][l])].PQ_Jpi_index==k)
								{
									pair_k_num++;
								}
							}
							block_mat[norp][pi][i].PP_in_basis_unnor[k][pack_k(p, q)] = block_mat[norp][pi][i].ove[pack_k(p,q)]*pair_k_num/2.0;
						}
					}
					else
					{
						for (int k = 0; k < N[norp] / 2; k++)
						{
							if (basis_block[norp][pi][i].r[p][k] > 0)
							{
								sr[k] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[p][k]].J, P_M[norp][basis_block[norp][pi][i].r[p][k]].pi, P_M[norp][basis_block[norp][pi][i].r[p][k]].M, P_M[norp][basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
							}
							else
							{
								sr[k] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[p][k]].J, P_M[norp][-basis_block[norp][pi][i].r[p][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[p][k]].M, P_M[norp][-basis_block[norp][pi][i].r[p][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[p][k]);
							}

							if (basis_block[norp][pi][i].r[q][k] > 0)
							{
								sr[k + N[norp] / 2+1] = PQ_mat_out(P_M[norp][basis_block[norp][pi][i].r[q][k]].J, P_M[norp][basis_block[norp][pi][i].r[q][k]].pi, P_M[norp][basis_block[norp][pi][i].r[q][k]].M, P_M[norp][basis_block[norp][pi][i].r[q][k]].PQ_Jpi_index,basis_block[norp][pi][i].r[q][k]);
							}
							else
							{
								sr[k + N[norp] / 2+1] = PQ_mat_out(P_M[norp][-basis_block[norp][pi][i].r[q][k]].J, P_M[norp][-basis_block[norp][pi][i].r[q][k]].pi, -P_M[norp][-basis_block[norp][pi][i].r[q][k]].M, P_M[norp][-basis_block[norp][pi][i].r[q][k]].PQ_Jpi_index, basis_block[norp][pi][i].r[q][k]);
							}
						}
						if(basis_block[norp][pi][i].r[p][N[norp]/2]>0)
						{
							int sp_index=basis_block[norp][pi][i].r[p][N[norp]/2]-1;
							sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,0);
						}
						else
						{
							int sp_index=-basis_block[norp][pi][i].r[p][N[norp]/2]-1;
							sr[N[norp]/2] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),0);
						}
						if(basis_block[norp][pi][i].r[q][N[norp]/2]>0)
						{
							int sp_index=basis_block[norp][pi][i].r[q][N[norp]/2]-1;
							sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, INT_MAX/4+sp_index+1,1);
						}
						else
						{
							int sp_index=-basis_block[norp][pi][i].r[q][N[norp]/2]-1;
							sr[(N[norp]/2)*2+1] = PQ_mat_out(sp_nljm[norp][sp_index].j, sp_nljm[norp][sp_index].l%2, -sp_nljm[norp][sp_index].m, sp_nljm[norp][sp_index].nlj, -INT_MAX/4-(sp_index+1),1);
						}
						sr_o_sort(sr,N[norp]/2);

						auto Q_Jp = PQ_mat_out(2, 0, 2, -1, J_pm_index[norp]);
						auto Q_Jm = PQ_mat_out(2, 0, -2, -1, -J_pm_index[norp]);
						block_mat[norp][pi][i].JpJm_unnor[pack_k(p, q)] = qq_o_cal(sr, Q_Jp, Q_Jm, N[norp] / 2, norp);
						
						PQ_mat_free(Q_Jp);
						PQ_mat_free(Q_Jm);

						for (int k = 0; k < sp_nlj[norp].size(); k++)
						{
							PQ_mat_struct Q0 = PQ_mat_out(0, 0, 0, -1, Q_ob_index[norp][k]);
							block_mat[norp][pi][i].spe_unnor[k][pack_k(p, q)] = q_o_cal(sr, Q0, N[norp] / 2, norp);
							PQ_mat_free(Q0);
						}
						for (int k = 0; k < PP_num[norp]; k++)
						{
							int J = PP_ham[norp][k].J;
							int pi_PP = PP_ham[norp][k].pi;
							int M_num = PP_ham[norp][k].PQ_index.size();
							for (int l = 0; l < M_num; l++)
							{
								int M = P_M[norp][PP_ham[norp][k].PQ_index[l]].M;
								int PQ_Jpi_index = P_M[norp][PP_ham[norp][k].PQ_index[l]].PQ_Jpi_index;
								auto P = PQ_mat_out(J, pi_PP, M, PQ_Jpi_index, PP_ham[norp][k].PQ_index[l]);
								double delta = pp_o_cal(sr, P, P, N[norp] / 2, norp);
								block_mat[norp][pi][i].PP_unnor[k][pack_k(p, q)] += delta;
								// block_mat_temp.PP[k][p+q*dim]+=pp_cal(sr,P,P,N[norp]/2,norp);
								PQ_mat_free(P);
								if (M != 0)
								{
									P = PQ_mat_out(J, pi_PP, -M, PQ_Jpi_index, -PP_ham[norp][k].PQ_index[l]);
									double delta = pp_o_cal(sr, P, P, N[norp] / 2, norp);

									block_mat[norp][pi][i].PP_unnor[k][pack_k(p, q)] += delta;
									// block_mat_temp.PP[k][p+q*dim]+=pp_cal(sr,P,P,N[norp]/2,norp);
									PQ_mat_free(P);
								}
							}
						}
						for (int k = 0; k < QQ_num[norp]; k++)
						{
							int J = QQ_ham[norp][k].J;
							int pi_QQ = QQ_ham[norp][k].pi;
							int M_num = QQ_ham[norp][k].PQ_index.size();
							for (int l = 0; l < M_num; l++)
							{
								int M = Q_M[norp][QQ_ham[norp][k].PQ_index[l] - INT_MAX / 2].M;
								int PQ_Jpi_index = Q_M[norp][QQ_ham[norp][k].PQ_index[l] - INT_MAX / 2].PQ_Jpi_index;
								auto Q = PQ_mat_out(J, pi_QQ, M, PQ_Jpi_index, QQ_ham[norp][k].PQ_index[l]);
								auto QT = PQ_mat_out(J, pi_QQ, -M, PQ_Jpi_index, -QQ_ham[norp][k].PQ_index[l]);
								double delta = qq_o_cal(sr, Q, QT, N[norp] / 2, norp);
								block_mat[norp][pi][i].QQ_unnor[k][pack_k(p, q)] += delta;
								// block_mat_temp.QQ[k][p+q*dim]+=qq_cal(sr,Q,QT,N[norp]/2,norp);
								if (M != 0)
								{
									delta = qq_o_cal(sr, QT, Q, N[norp] / 2, norp);
									block_mat[norp][pi][i].QQ_unnor[k][pack_k(p, q)] += delta;
									// block_mat_temp.QQ[k][p+q*dim]+=qq_cal(sr,QT,Q,N[norp]/2,norp);
								}
								PQ_mat_free(Q);
								PQ_mat_free(QT);
							}
						}
						for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
						{
							int pair_k_num=0;
							for(int l=0;l<N[norp]/2;l++)
							{
								if(P_M[norp][abs(basis_block[norp][pi][i].r[p][l])].PQ_Jpi_index==k)
								{
									pair_k_num++;
								}
								if(P_M[norp][abs(basis_block[norp][pi][i].r[q][l])].PQ_Jpi_index==k)
								{
									pair_k_num++;
								}
							}
							block_mat[norp][pi][i].PP_in_basis_unnor[k][pack_k(p, q)] = block_mat[norp][pi][i].ove[pack_k(p,q)]*pair_k_num/2.0;
						}
						for (int k = 0; k < sp_nlj[norp].size(); k++)
						{
							int C_num=0;
							if(sp_nljm[norp][abs(basis_block[norp][pi][i].r[p][N[norp]/2])-1].nlj==k)
							{
								C_num++;
							}
							if(sp_nljm[norp][abs(basis_block[norp][pi][i].r[q][N[norp]/2])-1].nlj==k)
							{
								C_num++;
							}
							block_mat[norp][pi][i].C_in_basis_unnor[k][pack_k(p, q)] = block_mat[norp][pi][i].ove[pack_k(p,q)]*C_num*0.5;
						}
					}
					
					
					cout<<block_mat[norp][pi][i].spe_unnor[1][pack_k(p, q)]<<' ';
				}
				cout<<endl;
			}
		}
		ham_mat_ort_nor(block_mat[norp][pi][i].JpJm,block_mat[norp][pi][i].JpJm_unnor, block_mat[norp][pi][i].vec, dim);
		for (int k = 0; k < sp_nlj[norp].size(); k++)
		{
			ham_mat_ort_nor(block_mat[norp][pi][i].spe[k],block_mat[norp][pi][i].spe_unnor[k], block_mat[norp][pi][i].vec, dim);
		}
		for (int k = 0; k < PP_num[norp]; k++)
		{
			ham_mat_ort_nor(block_mat[norp][pi][i].PP[k],block_mat[norp][pi][i].PP_unnor[k], block_mat[norp][pi][i].vec, dim);
		}
		for (int k = 0; k < QQ_num[norp]; k++)
		{
			ham_mat_ort_nor(block_mat[norp][pi][i].QQ[k],block_mat[norp][pi][i].QQ_unnor[k], block_mat[norp][pi][i].vec, dim);
		}
		for (int k = 0; k < P_Jpi_in_basis_num[norp]; k++)
		{
			ham_mat_ort_nor(block_mat[norp][pi][i].PP_in_basis[k],block_mat[norp][pi][i].PP_in_basis_unnor[k], block_mat[norp][pi][i].vec, dim);
		}
		if(N[norp]%2==1)
		{
			for (int k = 0; k < sp_nlj[norp].size(); k++)
			{
				ham_mat_ort_nor(block_mat[norp][pi][i].C_in_basis[k],block_mat[norp][pi][i].C_in_basis_unnor[k], block_mat[norp][pi][i].vec, dim);
			}
		}
		block_mat_ham_evalued[norp][pi][i]=true;
	}
}