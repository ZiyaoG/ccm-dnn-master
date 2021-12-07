% check the network trend
params = load('params.mat').params;
prd_dist = @(x) -uncer_func(x,params);

x = zeros([100,4]);
for i=1:100
    x(i,:) = [5 0.1*i 1  2    ];
end

d = prd_dist(x);
d1 = d(1,:);
d2 = d(2,:);



figure(1)
clf;
plot(d1)
hold on 
plot(d2)