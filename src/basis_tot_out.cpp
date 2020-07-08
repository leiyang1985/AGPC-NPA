struct basis_tot_block_struct
{
    int M[2];
    int pi[2];
    int block[2];
    int begin;
    int end;
    int dim;
};
vector<basis_tot_block_struct> basis_tot_block[2];
vector<array<int,6>> basis_tot[2];
int dim_tot[2]={0,0};
bool basis_tot_cmp(array<int,6> T1,array<int,6> T2)
{
    if(
        T1[0]<T2[0]||(T1[0]==T2[0]&&T1[3]<T2[3])
        ||
        (T1[0]==T2[0]&&T1[3]==T2[3]&&T1[1]<T2[1])
        ||
        (T1[0]==T2[0]&&T1[3]==T2[3]&&T1[1]==T2[1]&&T1[4]<T2[4])
        ||
        (T1[0]==T2[0]&&T1[3]==T2[3]&&T1[1]==T2[1]&&T1[4]==T2[4]&&T1[2]<T2[2])
        ||
        (T1[0]==T2[0]&&T1[3]==T2[3]&&T1[1]==T2[1]&&T1[4]==T2[4]&&T1[2]==T2[2]&&T1[5]<T2[5])
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}
void basis_tot_out(int M_tot)
{
    for(int pi_0=0;pi_0<2;pi_0++)
    {
        for(int pi_1=0;pi_1<2;pi_1++)
        {
            for(int i=0;i<basis_block[0][pi_0].size();i++)
            {
                for(int j=0;j<basis_block[1][pi_1].size();j++)
                {
                    if(basis_block[0][pi_0][i].M+basis_block[1][pi_1][j].M==M_tot)
                    {
                        int pi=(pi_0+pi_1)%2;
                        dim_tot[pi]+=basis_block[0][pi_0][i].r.size()*basis_block[1][pi_1][j].r.size();
                        for(int p=0;p<basis_block[0][pi_0][i].r.size();p++)
                        {
                            for(int q=0;q<basis_block[1][pi_1][j].r.size();q++)
                            {
                                array<int,6> one_basis {pi_0,i,p,pi_1,j,q};
                                basis_tot[pi].emplace_back(one_basis);
                            }
                        }
                    }
                }
            }
        }
    }
    cout<<"dim="<<dim_tot[0]<<' '<<dim_tot[1]<<endl;
    for(int pi=0;pi<2;pi++)
    {
        cout<<"go into basis_sort"<<endl;
        sort(basis_tot[pi].begin(),basis_tot[pi].end(),basis_tot_cmp);
        cout<<"go out of basis_sort"<<endl;
        // for(int i=0;i<dim_tot[pi];i++)
        // {
        //     for(int k=0;k<6;k++)
        //     {
        //         cout<<basis_tot[pi][i][k]<<' ';
        //     }
        //     cout<<endl;
        // }
        int begin=0;
        int end=0;
        int dim;
        int dim2;
        if(basis_tot[pi].size()>1)
        {
            for(auto i=basis_tot[pi].begin();i!=basis_tot[pi].end()-1;i++)
            {
                if((*i)[0]!=(*(i+1))[0]||(*i)[3]!=(*(i+1))[3]||(*i)[1]!=(*(i+1))[1]||(*i)[4]!=(*(i+1))[4])
                {
                    end=distance(basis_tot[pi].begin(),i+1);
                    cout<<"from "<<begin<<" to "<<end<<endl;
                    basis_tot_block_struct block_temp;
                    block_temp.pi[0]=(*i)[0];
                    block_temp.pi[1]=(*i)[3];
                    block_temp.block[0]=(*i)[1];
                    block_temp.block[1]=(*i)[4];
                    block_temp.M[0]=basis_block[0][block_temp.pi[0]][block_temp.block[0]].M;
                    block_temp.M[1]=basis_block[1][block_temp.pi[1]][block_temp.block[1]].M;
                    block_temp.begin=begin;
                    block_temp.end=end;
                    block_temp.dim=end-begin;
                    dim=end-begin;
                    dim2=dim*dim;
                    basis_tot_block[pi].emplace_back(block_temp);
                    begin=distance(basis_tot[pi].begin(),i+1);
                }
            }
        }
        if(basis_tot[pi].size()>0)
        {
            end=basis_tot[pi].size();
            cout<<"from "<<begin<<" to "<<end<<endl;
            basis_tot_block_struct block_temp;
            block_temp.pi[0]=basis_tot[pi][basis_tot[pi].size()-1][0];
            block_temp.pi[1]=basis_tot[pi][basis_tot[pi].size()-1][3];
            block_temp.block[0]=basis_tot[pi][basis_tot[pi].size()-1][1];
            block_temp.block[1]=basis_tot[pi][basis_tot[pi].size()-1][4];
            block_temp.M[0]=basis_block[0][block_temp.pi[0]][block_temp.block[0]].M;
            block_temp.M[1]=basis_block[1][block_temp.pi[1]][block_temp.block[1]].M;
            block_temp.begin=begin;
            block_temp.end=end;
            block_temp.dim=end-begin;
            dim=end-begin;
            dim2=dim*dim;
            basis_tot_block[pi].emplace_back(block_temp);
        }
    }
}
