% Solving the 2-D Poisson equation by the Finite Difference Method 
% Numerical scheme used is a second order central difference in space (5-point difference)
% https://es.mathworks.com/matlabcentral/fileexchange/38090-2d-poisson-equation

function p = poisson_solver_2D(nx,ny,niter,b,x,y,dx,dy)
%%
pn=zeros(nx,ny);                 %Preallocating p
% Initial Conditions
p=zeros(nx,ny);                  %Preallocating p
%Boundary conditions
p(:,1)=0;
p(:,ny)=0;
p(1,:)=0;                  
p(nx,:)=0;
%% solver
i=2:nx-1;
j=2:ny-1;
delta = 1e-6;
%Explicit iterative scheme with C.D in space (5-point difference)
for it=1:niter
    pn=p;
    p(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1)))-(b(i,j)*dx^2*dy^2))/(2*(dx^2+dy^2));
    % Ensure Boundary conditions 
    p(:,1)=0;
    p(:,ny)=0;
    p(1,:)=0;                  
    p(nx,:)=0;
    res = max(max(abs(p-pn)));
%     fprintf('iteration: %d; residual: %e \n',it,res);
    if res < delta 
        break
    end
end

% % Plotting the solution
% pmax = max(max(p));
% figure 
% h=surf(x,y,(p'),'EdgeColor','none');       
% shading interp
% % axis([-0.5 2.5 -0.5 2.5 -100 100])
% title({'2-D Poisson equation';['{\itNumber of iterations} = ',num2str(it)]})
% xlabel('Spatial co-ordinate (x) \rightarrow')
% ylabel('{\leftarrow} Spatial co-ordinate (y)')
% zlabel('Solution profile (P) \rightarrow')

end