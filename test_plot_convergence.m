%% test_plot_convergence

figure(10)
subplot(2,1,1)
plot([1:1:10],abs(it_track))
hold on
axis([1 10 0 0.05])

subplot(2,1,2)
semilogy([1:1:10],co2_track)
hold on
axis([1 10 0 220])