void tran_cal()
{
	double delta;

	for (int i = 0; i < tran.size(); i++)
	{
		int norp = tran[i].Q_norp;
		if (N[norp] <= 0)
		{
			continue;
		}
		int k = tran[i].Q_i;
		int Q_index = Q_M_tran_index[norp][k] - INT_MAX / 2;
		int her = Q_tran_her[norp][k];
		int pi = Q_M[norp][Q_index].pi;
		int kappa = Q_M[norp][Q_index].J;
		int pi_r = tran[i].pi_r;
		int pi_l = tran[i].pi_l;
		if ((pi_r + pi_l + pi) % 2 != 0)
		{
			continue;
		}

		int dim_l = dim[pi_l];
		int J_l = tran[i].J_l;
		int order_l = tran[i].i_l;
		int k_l;
		bool exist_l=false;
		for (int ii = 0; ii < state_num[pi_l]; ii++)
		{
			if (J[pi_l][ii] == J_l && order_l == order[pi_l][ii])
			{
				k_l = ii;
				exist_l=true;
				break;
			}
		}
		if(!exist_l)
		{
			continue;
		}
		double E_l = e[pi_l][k_l];

		int dim_r = dim[pi_r];
		int J_r = tran[i].J_r;
		int order_r = tran[i].i_r;
		int k_r;
		bool exist_r=false;
		for (int ii = 0; ii < state_num[pi_r]; ii++)
		{
			if (J[pi_r][ii] == J_r && order_r == order[pi_r][ii])
			{
				k_r = ii;
				exist_r=true;
				break;
			}
		}
		if(!exist_r)
		{
			continue;
		}
		double E_r = e[pi_r][k_r];
		if (fabs(cg_out(J_r, M_tot, kappa, 0, J_l, M_tot)) < 1e-12)
		{
			continue;
		}

		double res = 0;
		for(auto ii=basis_tot_block[pi_l].begin();ii!=basis_tot_block[pi_l].end();ii++)
		{
			int pi0_l=ii->pi[0];
			int pi1_l=ii->pi[1];
			int i0_l=ii->block[0];
			int j1_l=ii->block[1];
			int dim0_l=basis_block[0][pi0_l][i0_l].r.size();
			int dim1_l=basis_block[1][pi1_l][j1_l].r.size();
			int basis_tot_block_begin_l=ii->begin;
			int basis_tot_block_dim_l=ii->dim;
			int basis_tot_block_end_l=ii->end;
			for(auto jj=basis_tot_block[pi_r].begin();jj!=basis_tot_block[pi_r].end();jj++)
			{
				int pi0_r=jj->pi[0];
				int pi1_r=jj->pi[1];
				int i0_r=jj->block[0];
				int j1_r=jj->block[1];
				int dim0_r=basis_block[0][pi0_r][i0_r].r.size();
				int dim1_r=basis_block[1][pi1_r][j1_r].r.size();
				int basis_tot_block_begin_r=jj->begin;
				int basis_tot_block_dim_r=jj->dim;
				int basis_tot_block_end_r=jj->end;
				if (j1_l!=j1_r&&norp==0)
				{
					continue;
				}
				if (i0_l!=i0_r&&norp==1)
				{
					continue;
				}
				for(auto p=basis_tot[pi_l].begin()+basis_tot_block_begin_l;p!=basis_tot[pi_l].begin()+basis_tot_block_end_l;p++)
				{
					int p0_l=(*p)[2];
					int q1_l=(*p)[5];
					for(auto q=basis_tot[pi_r].begin()+basis_tot_block_begin_r;q!=basis_tot[pi_r].begin()+basis_tot_block_end_r;q++)
					{
						int p0_r=(*q)[2];
						int q1_r=(*q)[5];
						if(norp==0&&q1_l!=q1_r)
						{
							continue;
						}
						if(norp==1&&p0_l!=p0_r)
						{
							continue;
						}
						if(norp==0)
						{
							res+=v[pi_l][k_l][distance(basis_tot[pi_l].begin(),p)]*v[pi_r][k_r][distance(basis_tot[pi_r].begin(),q)]*tran_q_block[norp][pi0_l][pi0_r][i0_l][i0_r].mat[k][p0_l+p0_r*dim0_l];
						}
						if(norp==1)
						{
							res+=v[pi_l][k_l][distance(basis_tot[pi_l].begin(),p)]*v[pi_r][k_r][distance(basis_tot[pi_r].begin(),q)]*tran_q_block[norp][pi1_l][pi1_r][j1_l][j1_r].mat[k][q1_l+q1_r*dim1_l];
						}
					}
				}
			}
		}
		if (J_l == J_r && pi_l == pi_r && order_l == order_r)
		{
			cout << "<" << J_l << "_" << order_l << "^" << pi_l << ",M=" << J_l << ",e=" << E_l << "|Q_" << k << "^" << norp << "|>: <||~||>= " << res / cg_cal(J_l, M_tot, Q_M[norp][Q_index].J, 0, J_l, M_tot) << " CG*<||~||>= " << res / cg_cal(J_l, M_tot, Q_M[norp][Q_index].J, 0, J_l, M_tot) * cg_cal(J_l, J_l, Q_M[norp][Q_index].J, 0, J_l, J_l) << endl;
		}
		else
		{
			cout << "<" << J_l << "_" << order_l << "^" << pi_l << ",e_l=" << E_l << "||Q_" << k << "^" << norp << "||" << J_r << "_" << order_r << "^" << pi_r << ",e_r=" << E_r << ">: <||~||>= " << res / cg_cal(J_r, M_tot, Q_M[norp][Q_index].J, 0, J_l, M_tot) << " sqrt(hat J_l/ hat J_r)*<||~||>= " << res / cg_cal(J_r, M_tot, Q_M[norp][Q_index].J, 0, J_l, M_tot) * sqrt((J_l + 1) / double(J_r + 1)) << " D E= " << E_r - E_l << endl;
		}
	}
}