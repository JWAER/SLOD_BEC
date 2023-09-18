function ρvᵢ(U,W)
	Rho = zeros(size(U));
	
	Uconj = U';
	
@inbounds	for i = 1:length(W.Iptr)-1
		for idx = W.Iptr[i]:W.Iptr[i+1]-1;
			j = W.J[idx];
			k = W.K[idx];
			
			if(i==k); Rho[i]+=W.Val[idx]*real(U[j]*Uconj[k]); 
			elseif(j==k); Rho[i]+=W.Val[idx]*real(U[j]*Uconj[k]); Rho[j]+= 2*W.Val[idx]*real(U[i]*Uconj[k]);  
			elseif(i==j); Rho[i]+=2*W.Val[idx]*real(U[j]*Uconj[k]); Rho[k]+=W.Val[idx]*real(U[i]*Uconj[j]);   
			else
		
				
				Rho[i] += 2*W.Val[idx]*real(U[j]*Uconj[k]);
				Rho[j] +=  2*W.Val[idx]*real(U[i]*Uconj[k]);
				Rho[k] += 2*W.Val[idx]*real(U[i]*Uconj[j]);
				
			 end
			
			
		end
	end

	
	return Rho

end





function u²uvᵢ(W,U,Rho); #Rho here is P(u^2)

NL = similar(U); fill!(NL,0);

@inbounds	for i = 1:length(W.Iptr)-1
		for idx = W.Iptr[i]:W.Iptr[i+1]-1;
			j = W.J[idx];
			k = W.K[idx];
			
			if(i==k) NL[i]+= W.Val[idx]*(Rho[j]*U[k]);  #Symetric 
			elseif(j==k); NL[i]+=W.Val[idx]*(Rho[j]*U[k]);NL[j]+=W.Val[idx]*(Rho[i]*U[k]+Rho[k]*U[i]); #ok
			elseif(i==j); NL[i]+=W.Val[idx]*(Rho[j]*U[k]+Rho[k]*U[j]);NL[k] += W.Val[idx]*(Rho[i]*U[j]);
			else
				NL[i] += W.Val[idx]*(Rho[j]*U[k]+Rho[k]*U[j]);
				NL[j] += W.Val[idx]*(Rho[i]*U[k]+Rho[k]*U[i]);
				NL[k] += W.Val[idx]*(Rho[i]*U[j]+Rho[j]*U[i]);
			end
		end
	end

	return NL


end





