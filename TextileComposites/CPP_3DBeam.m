x1bar = [1,5.25e-5,-6.157e-22]
x2bar = [2,0.000168,-7.88e-22]
x1 = [1,0.0009475,-1.2e-21]
x2 = [2,0.000832,-2.62e-22]
b = x2 + x1
t = x2 - x1
bbar = x2bar + x1bar
tbar = x2bar - x1bar
exp = (tbar*dot(tbar, t) - t*dot(tbar,tbar))/(dot(tbar,tbar)*dot(t,t) - dot(tbar,t)^2)
dot(-(bbar-b),exp)
exp2 = (t*dot(tbar, t) - tbar*dot(t,t))/(dot(tbar,tbar)*dot(t,t) - dot(tbar,t)^2)
dot(bbar - b,exp2)