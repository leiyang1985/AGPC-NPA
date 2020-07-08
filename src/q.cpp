double q_cal(PQ_mat_struct* sr,PQ_mat_struct &Q,int N_ove,int norp)
{
	PQ_mat_struct sr_low[(N_ove)*2];
	node<PQ_mat_struct*,double> *ove_tree=NULL;
	for (int i = 0; i < N_ove; i++)
	//int i=0;
	{
		int coe_pq;
		sr_low[i]=PQ_con_cal(sr[i], Q,coe_pq);
		PQ_mat_cpy(sr_low, sr,  i);
		PQ_mat_cpy(sr_low + i+1, sr + i + 1,  (N_ove ) * 2 - i-1);
		double coe=coe_pq;
		sr_sort(sr_low,N_ove);
		ove_tree_search_add(ove_tree,sr_low,coe,sr_cmp,N_ove);
	}
	return tree_ove_op(ove_tree,N_ove,norp);
}

double q_o_cal(PQ_mat_struct* sr,PQ_mat_struct &Q,int N_ove,int norp)
{
	PQ_mat_struct sr_low[(N_ove)*2+2];
	node<PQ_mat_struct*,double> *ove_tree=NULL;
	for (int i = 0; i < N_ove; i++)
	//int i=0;
	{
		int coe_pq;
		sr_low[i]=PQ_con_cal(sr[i], Q,coe_pq);
		PQ_mat_cpy(sr_low, sr,  i);
		PQ_mat_cpy(sr_low + i+1, sr + i + 1,  (N_ove ) * 2 - i+1);
		double coe=coe_pq;
		sr_o_sort(sr_low,N_ove);
		ove_o_tree_search_add(ove_tree,sr_low,coe,sr_o_cmp,N_ove);
	}
	int coe_pq;
	sr_low[N_ove]=CQ_con_cal(sr[N_ove],Q,coe_pq);
	PQ_mat_cpy(sr_low, sr,  N_ove);
	PQ_mat_cpy(sr_low + N_ove+1, sr + N_ove + 1,  N_ove +1);
	double coe=coe_pq;
	sr_o_sort(sr_low,N_ove);
	ove_o_tree_search_add(ove_tree,sr_low,coe,sr_o_cmp,N_ove);
	return tree_ove_o_op(ove_tree,N_ove,norp);
}
