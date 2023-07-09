argList=argv();
numMat=str2num(argList{1})
numMotion=str2num(argList{2})

setXTicks

set(gcf,'Visible','off');
axis([0,numMotion+1]);

hold on

for i = 1:numMotion
  for j=1:numMat
      fileName=sprintf("NodeDisp.out.%d.%d",j,i)
      b=load(fileName);
      maxDisp=max(abs(b(:,2)));
      a(j,1)=i;
      a(j,2)=maxDisp;
      clear b;
   end
   plot(a(:,1),a(:,2),'*')
   plot(a(1,1),a(1,2),'or')
end

ylabel('Relative Roof Displacement')
xlabel('Earthquake')

setXTicks

print -djpg Disp.jpg

hold off
clf
hold on


axis([0,numMotion+1]);
xlabel('Motion');
ylabel('Max Accel');
for i = 1:numMotion
  for j=1:numMat
      fileName=sprintf("NodeAccel.out.%d.%d",j,i)
      b=load(fileName);
      maxDisp=max(abs(b(:,2)));
      a(j,1)=i;
      a(j,2)=maxDisp;
      clear b;
   end
   plot(a(:,1),a(:,2),'*')
   plot(a(1,1),a(1,2),'or')
end


ylabel('Absolute Roof Accelerations')
xlabel('Earthquake')

setXTicks
   

print -djpg Accel.jpg


hold off
clf
hold on


axis([0,numMotion+1]);
for i = 1:numMotion
  for j=1:numMat
      fileName=sprintf("NodeDrift.out.%d.%d",j,i)
      b=load(fileName);
      maxDisp=max(abs(b(:)));
      a(j,1)=i;
      a(j,2)=maxDisp;
      clear b;
   end
   plot(a(:,1),100*a(:,2),'*')
   plot(a(1,1),100*a(1,2),'or')
end


ylabel('Interstory Drift (percent)')
xlabel('Earthquake')

setXTicks
   

print -djpg Drift.jpg


hold off
clf
hold on


axis([0,numMotion+1]);
for i = 1:numMotion
  for j=1:numMat
      fileName=sprintf("NodeReaction.out.%d.%d",j,i)
      b=load(fileName);
      maxDisp=max(abs(sum(b,2)));
      maxDisp
      a(j,1)=i;
      a(j,2)=maxDisp;
      clear b;
   end
   a
   plot(a(:,1),100*a(:,2),'*')
   plot(a(1,1),100*a(1,2),'or')
end


ylabel('Base Shear')
xlabel('Earthquake')

setXTicks
   

print -djpg Reaction.jpg

quit
