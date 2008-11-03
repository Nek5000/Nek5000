function y = Efun(x)
	global lvl data ns
	y = x - amg_apply(data,lvl,ns,data.A{lvl}*x);
end
