function re=epsilon(k,param)

re=k.^2/(2*param.m)-param.mu;
% re=2*param.t*(1-cos(param.a*k))-param.mu;
end