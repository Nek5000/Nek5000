function [C,F]=coarse_par(A,tol)
	n = length(A(:,1));
	fprintf(1,'Computing S = I - D^{-1/2} A D^{-1/2} ...');
	D = diag(sparse(1./sqrt(diag(A))));
	S = abs(D*A*D);
	S = S - diag(diag(S));
	fprintf(1,' done, nnz(S) = %d.\n',nnz(S));
	vc = sparse(0*S(:,1));
	id = [1:n]';
	fprintf(1,'Running coarsening heuristic ...\n');
	while 1
		nc = ones(n,1)'*vc;
		vf = 1 - vc;
		g1 = vf.*(S*(vf.*(S*vf)));
		g2 = vf.*(S*(vf.*(S*g1)));
		g = vf.*sqrt(g2./(g1+vc));
		%pause; plot(g1); pause; plot(g2); pause;
		%plot(sqrt(g2 ./ (g1+vc))); pause;
		gd = vf .* (S*vf);
		[mg,mi] = max(g);
		fprintf(1,'  ratio = %g, n = %d, max Gershgorin radius = %g\n', ...
		  nc/n, nc, mg);
	  if(mg<=tol)
			if nc==0; vc(mi)=1; end
			break;
		end
		%[mg,mi] = max(gd);
		vc(mi)=1;
		vplot(vc); pause;
		%mask = g>tol;
		%ga = mask .* (g + .5*(S*(S*vc)));
		%ga = mask .* gd;
		%m = mxv_max(S,ga);
		%mask = mask & (ga-m>=0);
		%m = mxv_max(S,mask.*id);
		%mask = mask & (id-m>0);
    %vc = vc + mask;
	end
	C = find(vc);
	F = setdiff(id,C);
end
