function tot = wass_smooth(gamma,x,lambda,mu,sigma,m,lim)

tot = 0;
for i=1:m
q = integral(@(y)w_star(gamma,lambda(:,i),x,y,mu(i),sigma(i)),-lim,lim);
tot = tot + 1/m*(q);
end