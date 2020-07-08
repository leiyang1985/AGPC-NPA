double qq_cal(PQ_mat_struct* sr,PQ_mat_struct &Q1,PQ_mat_struct &Q2,int N_ove,int norp)
//the user needs to keep the sr serial workable for overlap calculation (M and parity conservation) before the calling of the overlap function
{
	PQ_mat_struct sr_low[N_ove*2];
	node<PQ_mat_struct*,double> *ove_tree=NULL;
	for (int i = 0; i < N_ove; i++)
	{
		for(int k=0;k<N_ove;k++)
		{
			if(i==k)
			{
				int coe_q1,coe_q2;
				PQ_mat_struct temp=PQ_con_cal(sr[i], Q1,coe_q1);
				sr_low[i]=PQ_con_cal(temp, Q2,coe_q2);
				PQ_mat_free(temp);
				PQ_mat_cpy(sr_low,sr,i);
				PQ_mat_cpy(sr_low+i+1,sr+i+1,2*N_ove-i-1);
				double coe=coe_q1*coe_q2;
				sr_sort(sr_low,N_ove);
				ove_tree_search_add(ove_tree,sr_low,coe,sr_cmp,N_ove);
			}
			else
			{
				int coe_q1,coe_q2;
				sr_low[i]=PQ_con_cal(sr[i], Q1,coe_q1);
				sr_low[k]=PQ_con_cal(sr[k], Q2,coe_q2);
				if(i<k)
				{
					PQ_mat_cpy(sr_low,sr,i);
					PQ_mat_cpy(sr_low+i+1,sr+i+1,k-i-1);
					PQ_mat_cpy(sr_low+k+1,sr+k+1,2*N_ove-k-1);
				}
				else
				{
					PQ_mat_cpy(sr_low,sr,k);
					PQ_mat_cpy(sr_low+k+1,sr+k+1,i-k-1);
					PQ_mat_cpy(sr_low+i+1,sr+i+1,2*N_ove-i-1);
				}
				double coe=coe_q1*coe_q2;
				sr_sort(sr_low,N_ove);
				ove_tree_search_add(ove_tree,sr_low,coe,sr_cmp,N_ove);
			}
		}
	}
	return tree_ove_op(ove_tree,N_ove,norp);
	//return 0;
}
double qq_o_cal(PQ_mat_struct* sr,PQ_mat_struct &Q1,PQ_mat_struct &Q2,int N_ove,int norp)
//the user needs to keep the sr serial workable for overlap calculation (M and parity conservation) before the calling of the overlap function
{
	PQ_mat_struct sr_low[N_ove*2+2];
	node<PQ_mat_struct*,double> *ove_tree=NULL;
	for (int i = 0; i < N_ove+1; i++)
	{
		for(int k=0;k<N_ove+1;k++)
		{
			if(i==k)
			{
				int coe_q1,coe_q2;
				if(i!=N_ove)
				{
					PQ_mat_struct temp=PQ_con_cal(sr[i], Q1,coe_q1);
					sr_low[i]=PQ_con_cal(temp, Q2,coe_q2);
					PQ_mat_free(temp);
				}
				else
				{
					PQ_mat_struct temp=CQ_con_cal(sr[i], Q1,coe_q1);
					sr_low[i]=CQ_con_cal(temp, Q2,coe_q2);
					PQ_mat_free(temp);
				}
				PQ_mat_cpy(sr_low,sr,i);
				PQ_mat_cpy(sr_low+i+1,sr+i+1,2*N_ove-i+1);
				double coe=coe_q1*coe_q2;
				sr_o_sort(sr_low,N_ove);
				ove_o_tree_search_add(ove_tree,sr_low,coe,sr_o_cmp,N_ove);
			}
			else
			{
				int coe_q1,coe_q2;
				if(i==N_ove)
				{
					sr_low[i]=CQ_con_cal(sr[i], Q1,coe_q1);
				}
				else
				{
					sr_low[i]=PQ_con_cal(sr[i], Q1,coe_q1);
				}
				
				if(k==N_ove)
				{
					sr_low[k]=CQ_con_cal(sr[k], Q2,coe_q2);
				}
				else
				{
					sr_low[k]=PQ_con_cal(sr[k], Q2,coe_q2);
				}
				if(i<k)
				{
					PQ_mat_cpy(sr_low,sr,i);
					PQ_mat_cpy(sr_low+i+1,sr+i+1,k-i-1);
					PQ_mat_cpy(sr_low+k+1,sr+k+1,2*N_ove-k+1);
				}
				else
				{
					PQ_mat_cpy(sr_low,sr,k);
					PQ_mat_cpy(sr_low+k+1,sr+k+1,i-k-1);
					PQ_mat_cpy(sr_low+i+1,sr+i+1,2*N_ove-i+1);
				}
				double coe=coe_q1*coe_q2;
				sr_o_sort(sr_low,N_ove);
				ove_o_tree_search_add(ove_tree,sr_low,coe,sr_o_cmp,N_ove);
			}
		}
	}
	return tree_ove_o_op(ove_tree,N_ove,norp);
	//return 0;
}
