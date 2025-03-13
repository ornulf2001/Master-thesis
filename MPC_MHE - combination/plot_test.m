x = 0:pi/100:2*pi;
y = sin(x);
figure(1)
plot(x,y)
hold on
plot(1)
hold off
text(3, 0.7, { 'Test 2'})