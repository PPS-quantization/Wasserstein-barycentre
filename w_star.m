function z = w_star(gamma,lambda,x,y,mu,sigma)

ci= pdist2(y',x','Euclidean').^2';
qy = normpdf(y,mu,sigma);

qy(qy < 1e-100) = 1e-100;

cc = exp(bsxfun(@plus, - (1/6.5)*ci, lambda)/gamma);
z = gamma*(log(1./qy)+ log(sum(cc))).*qy;

