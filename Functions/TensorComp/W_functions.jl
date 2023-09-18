

function getval_W(W,i,j,k)

	start = W.Iptr[i]
	finish = W.Iptr[i+1]-1;

	sub_j=searchsorted(W.J[start:finish], j)
	sub_k=searchsorted(W.K[start.+sub_j.-1], k)
	if(!isempty(sub_j)) return W.Val[(start+sub_j[1]-2).+sub_k]; end
	return 0.0

end



function addval_W(W,i,j,k,val)

	start = W.Iptr[i]
	finish = W.Iptr[i+1]-1;

	sub_j=searchsorted(W.J[start:finish], j)
	sub_k=searchsorted(W.K[start.+sub_j.-1], k)
	
	if(!isempty(sub_k));W.Val[(start+sub_j[1]-2).+sub_k].+=val; end
	return nothing
end


function print_idx_W(W,i,j,k)

	start = W.Iptr[i]
	finish = W.Iptr[i+1]-1;

	sub_j=searchsorted(W.J[start:finish], j)
	sub_k=searchsorted(W.K[start.+sub_j.-1], k)
	
	return (start+sub_j[1]-2).+sub_k
end




























