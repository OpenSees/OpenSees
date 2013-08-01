argList=argv();
load Node.out
plot(Node(:,2),Node(:,1));

ylabel('Force')
xlabel('Displacement')

print -djpg Disp.jpg

quit
