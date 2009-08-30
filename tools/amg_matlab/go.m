tol = input('target convergence factor [0.5]: ');
if length(tol) == 0; tol = 0.5; end
fprintf(1,'tol = %g\n',tol);

ctol = input('coarsening aggressiveness [0.9]: ');
if length(ctol) == 0; ctol = 0.9; end
fprintf(1,'ctol = %g\n',ctol);

[id A] = amg_import();
data = amg_setup(A, full(0*A(:,1)+1),ctol, tol);
Abottom = full(data.A{length(data.n)});
fprintf(1,'A_%d = %g (',length(data.n),Abottom);
if Abottom<1e-9
  nullspace = 1;
else
	nullspace = 0;
	fprintf(1,'no ');
end
fprintf(1,'nullspace)\n');
amg_export(data,id,nullspace);
