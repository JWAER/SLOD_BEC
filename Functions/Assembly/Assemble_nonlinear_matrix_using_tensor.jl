function ρφᵢφⱼ(Rho,W,Wprim,LTRI)


nzMLOD = length(LTRI.nzval);

N = length(Rho);
nzval = zeros(nzMLOD);
iit = 1;
col = zeros(N); #is actually row, but symmetry allows for it


	@inbounds for i = 1:length(W.Iptr)-1
		to_do = W.Iptr[i]:W.Iptr[i+1]-1;
		
		for idx = to_do;
			j = W.J[idx];
			k = W.K[idx];
			
			if(i==k) col[j] += W.Val[idx]*(Rho[j]);
			elseif(j==k); col[k] += W.Val[idx]*Rho[j]; 
			elseif(i==j); col[k] += W.Val[idx]*Rho[j]; col[j]+= W.Val[idx]*Rho[k];
			
			else
				col[k] += W.Val[idx]*Rho[j];
				col[j] += W.Val[idx]*Rho[k] 
			end
		end
		to_do = Wprim.Iptr[i]:Wprim.Iptr[i+1]-1;
		
		for idx = to_do
			jprim = Wprim.J[idx];
			kprim = Wprim.K[idx];
			col[jprim]+=Wprim.Val[idx]*Rho[kprim];
		end
		
		idx = LTRI.colptr[i]:LTRI.colptr[i+1]-1;
		
		idx_col_prealloc = LTRI.rowval[idx];
		nzval[idx] = col[idx_col_prealloc]; col[idx_col_prealloc].=0.#fill!(col,0);
		
		
	end

	@inbounds for c = 1:N
		for idx in nzrange(LTRI,c);
			if(c==LTRI.rowval[idx]) nzval[idx]/=2.;end
		end
	end
	
	
	M_NL_loc= SparseMatrixCSC(N,N,LTRI.colptr,LTRI.rowval,nzval)

	M_NL_loc += M_NL_loc';


	
	return M_NL_loc


end






















