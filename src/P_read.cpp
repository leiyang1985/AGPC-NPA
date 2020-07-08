vector<int> P_in_basis[2];
int P_Jpi_in_basis_num[2];
vector<int> P_num_max[2];//P_num_max=-1 means no limit of pair num
void P_read()
{
    FILE *fp;
	char filename[80]="P.dat";
	if((fp=fopen(filename,"r"))==NULL)
	{
		cout<<"open file error"<<endl;
	}
	else
	{
		int P_num_in_file;
		PQ_Jpi_struct P_Jpi_temp;
		PQ_M_struct P_mat_temp;
		int num_max_temp;
		for(int norp=0;norp<2;norp++)
		{
			fscanf(fp,"%d",&(P_num_in_file));
			for(int i=0;i<P_num_in_file;i++)
			{
				P_Jpi_temp.yq=new double [sp_nlj_dim2[norp]];
				fscanf(fp,"%d",&(P_Jpi_temp.J));
				fscanf(fp,"%d",&(P_Jpi_temp.pi));
				fscanf(fp,"%d",&num_max_temp);
				P_num_max[norp].emplace_back(num_max_temp);
				for(int k=0;k<sp_nlj_dim2[norp];k++)
				{
					fscanf(fp,"%lf", P_Jpi_temp.yq+k);
				}
				P_Jpi[norp].emplace_back(P_Jpi_temp);
				for(int M=0;M<=P_Jpi_temp.J;M+=2)
				{
					int P_M_index=PJpi_to_PM(P_Jpi[norp].size()-1,M,norp);
					P_in_basis[norp].emplace_back(P_M_index);
					if(M!=0)
					{
						P_in_basis[norp].emplace_back(-P_M_index);
					}
				}
			}
			P_Jpi_in_basis_num[norp]=P_Jpi[norp].size();
			//pow_gamma_out(norp);
		}
        fclose(fp);
    }
}
