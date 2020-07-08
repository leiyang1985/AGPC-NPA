int state_num[2];
int dim[2];
int *J[2];
int *order[2];
double *e[2];
double **v[2];
void vec_read()
{
    int temp_d;
    FILE *fp_vec=fopen("vec.bin","rb");
    fread(&M_tot,sizeof(int),1,fp_vec);
    for(int pi=0;pi<2;pi++)
    {
        fread(state_num+pi,sizeof(int),1,fp_vec);
        fread(dim+pi,sizeof(int),1,fp_vec);
        if(state_num<=0)
        {
            continue;
        }
        J[pi]=new int [state_num[pi]];
        
        order[pi]=new int [state_num[pi]];
        
        e[pi]=new double [state_num[pi]];
        
        v[pi]=new double *[state_num[pi]];
        for(int k=0;k<state_num[pi];k++)
        {
            fread(J[pi]+k,sizeof(int),1,fp_vec);
            fread(order[pi]+k,sizeof(int),1,fp_vec);
            fread(&temp_d,sizeof(int),1,fp_vec);
            fread(e[pi]+k,sizeof(double),1,fp_vec);
            v[pi][k]=new double [dim[pi]];
            fread(v[pi][k],sizeof(double),dim[pi],fp_vec);
        }
    }
    fclose(fp_vec);
}