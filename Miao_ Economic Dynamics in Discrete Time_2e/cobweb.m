
function cobweb(f,x0,N)
% generate the cobweb plot associated with
% the orbits x_n+1=f(x_n).
% N is the number of iterates, 
% x0 is the initial points.
% use @f to pass function ...



% turn hold on to gather up all plots in one
%hold on;
%plot(x,y,'k'); % plot the function
%plot(x,x,'r'); % plot the straight line
x(1)=x0; % plot orbit starting at x0
x(2)=f(x(1));
hold on
line([x(1),x(1)],[0,x(2)]);
line([x(1),x(2)],[x(2),x(2)])
deltax=x(2)-x(1);
quiver(x(1),x(2),deltax,0)
for i=2:N
    x(i+1)=f(x(i)); 
    line([x(i),x(i)],[x(i),x(i+1)]);
    quiver(x(i),x(i),0,x(i+1)-x(i));
    line([x(i),x(i+1)],[x(i+1),x(i+1)]);
    quiver(x(i),x(i+1),x(i+1)-x(i),0);
end

hold off;
