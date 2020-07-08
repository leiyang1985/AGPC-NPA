int N[2];
struct nlj_struct
{
	int j;
	int l;
	int N;
	bool unpair;
};
int sp_nlj_dim2[2];
vector<nlj_struct> sp_nlj[2];
struct nljm_struct
{
	int j;
	int l;
	int N;
	int nlj;
	int m;
};

int sp_nljm_dim2[2];
vector<nljm_struct> sp_nljm[2];
void sps_read()
{
    FILE *fp;
	char filename[80]="sps.dat";
	int temp_i;
	if((fp=fopen(filename,"r"))==NULL)
	{
		cout<<"open file error"<<endl;
	}
	else
	{
		for(int norp=0;norp<2;norp++)
		{
			int sp_nlj_dim;
			fscanf(fp,"%d",&(N[norp]));
			fscanf(fp,"%d",&(sp_nlj_dim));
			sp_nlj_dim2[norp]=sp_nlj_dim*sp_nlj_dim;
			nlj_struct temp_nlj[sp_nlj_dim];
			for(int i=0;i<sp_nlj_dim;i++)
			{
				fscanf(fp,"%d",&(temp_nlj[i].j));
			}
			for(int i=0;i<sp_nlj_dim;i++)
			{
				fscanf(fp,"%d",&(temp_nlj[i].l));
				temp_nlj[i].l/=2;
			}
			for(int i=0;i<sp_nlj_dim;i++)
			{
				fscanf(fp,"%d",&(temp_nlj[i].N));
			}
			for(int i=0;i<sp_nlj_dim;i++)
			{
				fscanf(fp,"%d",&temp_i);
				if(temp_i==1)
				{
					temp_nlj[i].unpair=true;
				}
				else
				{
					temp_nlj[i].unpair=false;
				}
			}
			sp_nlj[norp].insert(sp_nlj[norp].begin(),temp_nlj,temp_nlj+sp_nlj_dim);
			nljm_struct temp_nljm;
			for(vector<nlj_struct>::iterator i=sp_nlj[norp].begin();i!=sp_nlj[norp].end();i++)
			{
				for(temp_nljm.m=1;temp_nljm.m<=i->j;temp_nljm.m+=2)
				{
					temp_nljm.N=i->N;
					temp_nljm.l=i->l;
					temp_nljm.j=i->j;
					temp_nljm.nlj=distance(sp_nlj[norp].begin(),i);
					sp_nljm[norp].emplace(sp_nljm[norp].begin(),temp_nljm);
					temp_nljm.m=-temp_nljm.m;
					sp_nljm[norp].emplace_back(temp_nljm);
					temp_nljm.m=-temp_nljm.m;
				}
			}
			sp_nljm_dim2[norp]=sp_nljm[norp].size()*sp_nljm[norp].size();
		}
        fclose(fp);
    }
}