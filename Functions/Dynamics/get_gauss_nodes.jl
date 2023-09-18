function get_gauss_nodes(q,a,b)

if(q==1);
       tau0 = zeros(1);
        w0 = zeros(1);
        tau0[1] = 0;
        w0[1] = 2;
        

elseif(q==2);
        tau0 = [-sqrt(1/3), sqrt(1/3)];
        w0 = [1, 1];  
elseif(q==3)
  	tau0 = [-sqrt(3/5), 0,  sqrt(3/5)];
        w0 = [5/9, 8/9, 5/9]; 
elseif(q==4)
        tau0 = [-sqrt((3/7) + (2/7)*sqrt(6/5)), 
                -sqrt((3/7) - (2/7)*sqrt(6/5)), 
                sqrt((3/7) - (2/7)*sqrt(6/5)),
                sqrt((3/7) + (2/7)*sqrt(6/5)) ];
        w0 = [(18 - sqrt(30))/36,  
              (18 + sqrt(30))/36,  
              (18 + sqrt(30))/36, 
              (18 - sqrt(30))/36 ];  
elseif(q==6)
        tau0 = [ -0.932469514203152, 
                 -0.661209386466265, 
                 -0.238619186083197, 
                  0.238619186083197, 
                  0.661209386466265, 
                  0.932469514203152 ];
        w0 = [ 0.171324492379170, 
               0.360761573048139, 
               0.467913934572691, 
               0.467913934572691, 
               0.360761573048139, 
               0.171324492379170 ];  
        


end

	tau = tau0.*(b-a)/2 .+ (a+b)/2;
	w = w0.*(b-a)/2;

return tau,w
end


function lagrange_basis(q,τ)
if q == 2
li = [x-> (x-τ[2])./(τ[1]-τ[2]),  x-> (x-τ[1])./(τ[2] - τ[1]) ];


li_hat = [ x-> (x-τ[1]).*(x-τ[2])./(τ[1]*τ[2]); 
           x-> x.*(x-τ[2])./(τ[1]^2 - τ[1]*τ[2]); 
           x-> x.*(x-τ[1])./(τ[2]^2 - τ[1]*τ[2]) ];



# Derivatives of Lagrange polynomials
dli = [ x-> 1/(τ[1]-τ[2]); x-> 1/(τ[2] - τ[1]) ];


elseif( q== 3)

li = [ x-> (x-τ[i%q+1])*(x-τ[(i+1)%q+1])/(τ[i]-τ[i%q+1])/(τ[i]-τ[(i+1)%q+1]) for i = 1:q];

   
li_hat = [ x-> -((x-τ[1]).*(x-τ[2]).*(x-τ[3]))./(τ[1].*τ[2].*τ[3]); 
           x-> (x.*(x-τ[2]).*(x-τ[3]))./(τ[1].*(τ[1]-τ[2]).*(τ[1]-τ[3])); 
           x-> (x.*(x-τ[1]).*(x-τ[3]))./(τ[2].*(τ[2]-τ[1]).*(τ[2]-τ[3]));  
           x-> (x.*(x-τ[1]).*(x-τ[2]))./(τ[3].*(τ[3]-τ[1]).*(τ[3]-τ[2])) ];
            
dli = [ x-> -(τ[2]+τ[3]-2*x)./((τ[1]-τ[2])*(τ[1]-τ[3])); 
        x-> -(τ[1]+τ[3]-2*x)./((τ[2]-τ[1])*(τ[2]-τ[3])); 
        x-> -(τ[1]+τ[2]-2*x)./((τ[3]-τ[1])*(τ[3]-τ[2])) ];
  
end

return li, li_hat,dli


end