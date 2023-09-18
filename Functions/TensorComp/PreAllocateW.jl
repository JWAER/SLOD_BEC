
function PreAllocateW(M,prealloc=1);
SpaceDim = size(M,1);
W_Iptr = Array{Int64}(undef,SpaceDim+1)

W_Iptr[1] = 1;


N = size(M,1);

IT=0;

for basis = 1:N
        it = 0

        support = M[:,basis].nzind;
        n_loc = length(support)
        for j = 1:n_loc
                if(support[j]>=basis)
                        for k = j:n_loc
                        if(support[k]>=support[j])
                                        if(abs(M[support[j],support[k]])!=0.)
                                                it+=1; 
                                        end

                        end
                        end
                end
        end

        IT+=it;
        W_Iptr[basis+1] = IT+1; 

end


W_J = zeros(Int64,IT)
W_K = zeros(Int64,IT)
W_Val = zeros(IT);

it=0;


for basis = 1:N


        support = M[:,basis].nzind;
        n_loc = length(support)
        #println("n_loc ", n_loc)
        for j = 1:n_loc
                if(support[j]>=basis)
                        for k = j:n_loc
                        if(support[k]>=support[j])
                                        if(abs(M[support[j],support[k]])!=0.)
                                                it+=1; 
                                                W_J[it] = support[j]
                                                W_K[it] = support[k]
                                        end

                        end
                        end
                end
        end

        #println(basis)
#       fill!(mem,0);
end



return Tensor(W_Iptr,W_J,W_K,zeros(IT));
end



