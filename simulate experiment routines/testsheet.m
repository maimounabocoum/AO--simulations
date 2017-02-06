%% test sheet

w0 = 6e6;


dt = 1/(6e6);

t = (-1000:1000)*dt;

s = exp(-t.^2/(10e-6)^2).*cos(w0*t);


plot(t,s)