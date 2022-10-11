xx = randn(30,240);
tic;
gg1 = run_timeSeries2mat(xx,3,[0.0625,0.125],30,0,'wtc');
toc;
tic;
gg2 = run_timeSeries2matFastCohi(xx,3,[0.0625,0.125],30,0,'wtc');
toc;
diff = 0;
for i = 1:8
    diff = diff+sum(sum(abs(gg1{i}-gg2{i})));
end