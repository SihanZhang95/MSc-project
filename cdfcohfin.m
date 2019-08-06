function h=cdfcohfin(n, p, x, gam) %% % computes CDF of estimated mag sqd multiple coherence \gamh % Uses Gurland's approach of a finite sum of weighted incomplete betas% P(gamh < x) = \sum_{l=0}^{n-p+1} b_l I_z(p+l-1,n-p+1)%% n 		no. complex degrees of freedom% p			no. series% x			defines P(gamh <x)% gam		actual mag sqd multiple coherence%% h			prob%k=n-p+1; d2=k; % second dofd1=p-1;b0=(1-gam).^k;z=x.*(1-gam)./(1-gam.*x);t=b0.*betainc(z,d1,d2);h=t;for j=1:k	d1=p-1+j;% first dof	b1=((k-(j-1))./j).*b0;	b1=b1.*gam./(1-gam);	t=b1.*betainc(z,d1,d2);	b0=b1;	h=h+t;end	