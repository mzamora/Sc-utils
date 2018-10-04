function zb=get_wellmixed_zb(zi,qt,tl)

def=def_cnstnts; def.psrf=101000;

for it=1:length(zi) %if there's a time series of zi
    z=0:1:zi(it); 
    t_dry=nan(size(z)); t=t_dry; qsat=t_dry;
    p=def.psrf-1.2*def.g*z;
    pi=(p/def.p00).^def.rcp;

    %dry temperature
    t_dry=tl(it).*pi;
    qsat=rslf(p,t_dry)*1000;

    %zb is ay saturation of the dry temperature
    izb=find(qsat<qt(it),1);
    zb(it)=z(izb);

    %to get ql and T(z), extra steps are needed
end

end