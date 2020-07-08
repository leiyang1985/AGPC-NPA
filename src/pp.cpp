double pp_cal(PQ_mat_struct* sr,PQ_mat_struct &P1,PQ_mat_struct &P2,int N_ove,int norp)
//the user needs to keep the sr serial workable for overlap calculation (M and parity conservation) before the calling of the overlap function
{
	PQ_mat_struct sr_low[(N_ove)*2];
	node<PQ_mat_struct*,double> *ove_tree=NULL;
	for (int i = 0; i < N_ove; i++)
	//int i=0;
	{
		// double phi =1;
		double phi= PP_dot_cal(sr[i], P1, norp);

		if(fabs(phi)>1e-13)
		{
			PQ_mat_cpy(sr_low, sr,  i);
			PQ_mat_cpy(sr_low+i,&P2,1);
			PQ_mat_cpy(sr_low + i+1, sr + i + 1,  (N_ove ) * 2 - i-1);
			sr_sort(sr_low,N_ove);
			ove_tree_search_add(ove_tree,sr_low,phi,sr_cmp,N_ove);
		}		
	}
	int coe_pq;
	for (int i = 1; i < N_ove; i++)
	{
		auto Q = PP_con_cal(sr[i], P1);
		// cout<<"new "<<all<<" Q"<<endl;
		// PQ_mat_struct Q;
		// Q.val=new int [1];
		// Q.val[0]=0;
		// Q.ser=new int *[1];
		// Q.ser[0]=NULL;
		
		for (int k = 0; k < i; k++)
		{
			PQ_mat_cpy(sr_low, sr,  k);
			sr_low[k]=PQ_con_cal(sr[k],Q,coe_pq);
			PQ_mat_cpy(sr_low + k+1, sr + k + 1, i-k-1);
			PQ_mat_cpy(sr_low+i,&P2,1);
			PQ_mat_cpy(sr_low+i+1,sr+i+1,2*N_ove-i-1);
			sr_sort(sr_low,N_ove);
			double coe=coe_pq;
			ove_tree_search_add(ove_tree,sr_low,coe,sr_cmp,N_ove);
		}
		PQ_mat_free(Q);
	}
	return tree_ove_op(ove_tree,N_ove,norp);
	//return 0;
}
double pp_o_cal(PQ_mat_struct* sr,PQ_mat_struct &P1,PQ_mat_struct &P2,int N_ove,int norp)
//the user needs to keep the sr serial workable for overlap calculation (M and parity conservation) before the calling of the overlap function
{
	PQ_mat_struct sr_low[(N_ove)*2+2];
	node<PQ_mat_struct*,double> *ove_tree=NULL;
	for (int i = 0; i < N_ove; i++)
	//int i=0;
	{
		// double phi =1;
		double phi= PP_dot_cal(sr[i], P1, norp);

		if(fabs(phi)>1e-13)
		{
			PQ_mat_cpy(sr_low, sr,  i);
			PQ_mat_cpy(sr_low+i,&P2,1);
			PQ_mat_cpy(sr_low + i+1, sr + i + 1,  (N_ove ) * 2 - i+1);
			sr_o_sort(sr_low,N_ove);
			ove_o_tree_search_add(ove_tree,sr_low,phi,sr_o_cmp,N_ove);
		}		
	}
	int coe_pq;
	for (int i = 0; i < N_ove; i++)
	{
		auto Q = PP_con_cal(sr[i], P1);
		// cout<<"new "<<all<<" Q"<<endl;
		// PQ_mat_struct Q;
		// Q.val=new int [1];
		// Q.val[0]=0;
		// Q.ser=new int *[1];
		// Q.ser[0]=NULL;
		
		for (int k = 0; k < i; k++)
		{
			PQ_mat_cpy(sr_low, sr,  k);
			sr_low[k]=PQ_con_cal(sr[k],Q,coe_pq);
			PQ_mat_cpy(sr_low + k+1, sr + k + 1, i-k-1);
			PQ_mat_cpy(sr_low+i,&P2,1);
			PQ_mat_cpy(sr_low+i+1,sr+i+1,2*N_ove-i+1);
			sr_o_sort(sr_low,N_ove);
			double coe=coe_pq;
			ove_o_tree_search_add(ove_tree,sr_low,coe,sr_o_cmp,N_ove);
		}
		
		PQ_mat_cpy(sr_low, sr,  i);
		PQ_mat_cpy(sr_low+i,&P2,1);
		PQ_mat_cpy(sr_low+i+1,sr+i+1,N_ove-i-1);
		sr_low[N_ove]=CQ_con_cal(sr[N_ove],Q,coe_pq);
		PQ_mat_cpy(sr_low+N_ove+1,sr+N_ove+1,N_ove+1);
		sr_o_sort(sr_low,N_ove);
		double coe=coe_pq;
		ove_o_tree_search_add(ove_tree,sr_low,coe,sr_o_cmp,N_ove);
		PQ_mat_free(Q);
	}
	return tree_ove_o_op(ove_tree,N_ove,norp);
	//return 0;
}
