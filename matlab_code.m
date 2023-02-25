%%Perri Alessandro - s4476726
%%Assignment for Control of Linear Multivariable course

%First Request: Write the state equations of the system
%Let's define the needed matrices to properly represent our system:

A= [5 0 ; 2 3];
B=[1 0;0 2];
C=[0 1];
n=rank(A);
%now let's define the weight matrices:
V=[10 0; 0 10];
P=[1 0; 0 2];


%Second request: determine the stationary state-feedback control law 
%structure u(t) = -Lx(t)

[L,K, e]=lqr(A,B,V,P,0); %respectively L (optimal gain matrix), K (solution
                         % of Riccati equation), e (closed loop poles)

%Third Request: Check whether this control law ensures that the closed-loop
%system is asymptotically stable

%Let's check if the controllability matrix is full rank
ctrlMatrix=ctrb(A,B);
if(rank(ctrlMatrix)==n) %if the rank of the controllability matrix is equal
                        % to n then the system is fully controllable
    disp("The system is fully controllable")
end

H=chol(V);

Q=[H,H*A];

if(rank(Q)==n)
    disp("I've been able to find a matrix H such that V=transpose(H)*H and " + ...
        "such that it satisfies the condition rank(q)=n, so the system is asymptotically stable ")
end



%Fourth Request: Assume now that the state of the system is not fully accessible (but only the output y is). Design an 
% asymptotic observer of the state choosing the poles of the observer so that the observer dynamics is 
% considerably “faster” than the dynamics of the original closed-loop system (e.g., ensuring that the 
% time constants of the observer are an order of magnitude smaller than those of that system).

syms g1 g2
G=[g1;g2];
F=A-G*C;

obsMatrix=obsv(A,C); %observability matrix

if(rank(obsMatrix)==n) %if the rank of the observability matrix is equal to
                       % n then the system is fully observable
    disp("The system is fully observable")
end

%choose the values of the eigenvalues
lambda1=-10;
lambda2=-20;
 
eq1=det(lambda1*eye(2)-F)
eq2=det(lambda2*eye(2)-F)

[g1,g2]=solve(eq1==0,eq2==0) %find the values for g1 and g2


%Sixth Request: Determine the closed-loop transfer function 𝑻(𝒔) = 𝒀(𝒔)/𝑹(𝒔) (assuming a very fast convergence of the 
%observed system state to the true system state) and analyze the frequency response (for different frequencies 
%of the sinusoidal input).

%From the fifth point we know that u(t)=-L*x(t)+transpose([1 0])*r(t)

A_new=[A-B*L]
B_new=[B*[1;0]]

[b,a]=ss2tf(A_new,B_new,C,0)

sys= tf(b,a)

bode(sys)
grid on
[mag,phase]=bode(sys);

%Now we analyze the frequency response by considering different values of
%frequency:

%w=0.1, A=1
w=0.1
f=j*w
mod=abs(evalfr(sys,f))

%w=1, A=1
w=1
f=j*w
mod=abs(evalfr(sys,f))


%w=10, A=1
w=10
f=j*w
mod=abs(evalfr(sys,f))





