clear all;
close all;
clc;
disp('LIFT,DRAG,PRESSURE,TEMPERATURE,Cl and CD FOR DOUBLE WEDGE AIRFOIL FOR')
disp('SUPERSONIC FLOW ')
disp(' ')
disp('ENTER THE FOLLOWING DATA')
t=input('Enter the value of maximum thickness of double-wedge airfoil in meters: ');
c=input('Enter the value of chord length of double-wedge airfoil in meters: ');
M=input('Enter the free stream Mach number: ');
degrees=input('Enter the anti-clockwise angle of attack in degrees: ');
P=input('Enter the free stream static pressure in Pascal: ');
gamma=input('Enter the ratio of specific heats: ');
rho=input('Enter the free stream density in kg/m^3: ');
a=abs(degrees./57.2957795);
theta=atan(t./c);
cs=sqrt(gamma.*P./rho);
l=0.5.*(sqrt(t.^2+c.^2));
Dv=((M.^2-1).^2-3.*(1+((gamma-1)./2).*(M.^2)).*(1+((gamma+1)./2).*(M.^2)).*((tan(abs(theta)+abs(a))).^2));
if(Dv<0)
    disp(' ')
    disp(' ')
    disp('Detached shock exist or invalid data entry: Computation not possible !');
else
    Q2=abs(theta)-abs(a);
    disp(' ')
    disp('SOLUTION')
    disp('Valid data entered for free stream condition')
    if(Q2>0)
        disp('Fore upper region of airfoil undergoes a shock wave')
        f = @(b) tan(Q2)-(2.*cot(b))*((M.^2.*sin(b).^2-1)/(M.^2.*(gamma+cos(2.*b))+2));
        options = optimset('Display','off');
        beta2 = fsolve(f,0.001,options);
        P2 = P.*(1+((2.*gamma)./(gamma+1)).*((M.*sin(beta2)).^2-1));
        M2 = (sqrt(((M.*sin(beta2).^2)+(2./(gamma-1)))./(((2.*gamma)./(gamma-1)).*M.*sin(beta2)-1)))./(sin(beta2-Q2));
    else
        disp('Fore upper region of airfoil undergoes a expansion')
        syms IM
        g = @(m) abs(Q2)+sqrt((gamma+1)./(gamma-1)).*(atan(sqrt(((gamma-1)./(gamma+1)).*(M.^2-1)))) - ...
                 atan(sqrt(M.^2-1))-sqrt((gamma+1)./(gamma-1)).*(atan(sqrt(((gamma-1)./(gamma+1)).*(m.^2-1))))+atan(sqrt(m.^2-1));
        options = optimset('Display','off');
        M2 = fsolve(g,M,options);
        M2 = abs(M2);
        P2 = P.*((1+((gamma-1)./2).*(M.^2))./(1+((gamma-1)./2).*(M2.^2))).^(r./(r-1));
    end
    
    disp('Fore lower region of airfoil under goes a shock wave')
    Q3=abs(theta)+abs(a);
    f = @(b) tan(Q3)-(2.*cot(b))*((M.^2.*sin(b).^2-1)/(M.^2.*(gamma+cos(2.*b))+2));
    options = optimset('Display','off');
    beta3 = fsolve(f,0.001,options);
    P3 = P.*(1+((2.*gamma)./(gamma+1)).*((M.*sin(beta3)).^2-1));
    M3 = (sqrt(((M.*sin(beta3).^2)+(2./(gamma-1)))./(((2.*gamma)./(gamma-1)).*M.*sin(beta3)-1)))./(sin(beta3-Q3));
    
    disp('Aft upper region of aerofoil under goes expansion')
    Q4=2.*abs(theta);
    g = @(m) Q4 + sqrt((gamma+1)./(gamma-1)).*(atan(sqrt(((gamma-1)./(gamma+1)).*(M2.^2-1)))) - ...
             atan(sqrt(M2.^2-1))-sqrt((gamma+1)./(gamma-1)).*(atan(sqrt(((gamma-1)./(gamma+1)).*(m.^2-1))))+atan(sqrt(m.^2-1));
    M4 = fsolve(g,M2,options);
    M4 = abs(M4);
    P4 = P2.*((1+((gamma-1)./2).*(M2.^2))./(1+((gamma-1)./2).*(M4.^2))).^(gamma./(gamma-1));
    
    disp('Aft lower region of aerofoil under goes expansion')
    Q5=2.*abs(theta);
    g = @(m) Q5 + sqrt((gamma+1)./(gamma-1)).*(atan(sqrt(((gamma-1)./(gamma+1)).*(M3.^2-1)))) - ...
             atan(sqrt(M3.^2-1))-sqrt((gamma+1)./(gamma-1)).*(atan(sqrt(((gamma-1)./(gamma+1)).*(m.^2-1))))+atan(sqrt(m.^2-1));
    M5 = fsolve(g,M3,options);
    M5 = abs(M5);
    P5 = P3.*((1+((gamma-1)./2).*(M3.^2))./(1+((gamma-1)./2).*(M5.^2))).^(gamma./(gamma-1));
end

L=l.*((P5-P2).*cos(abs(theta)-abs(a))+(P3-P4).*cos(abs(theta)+abs(a)));
D=l.*((P5-P2).*sin(abs(theta)-abs(a))+(P3-P4).*sin(abs(theta)+abs(a)));
V=M.*cs;
dp=(0.5).*rho.*(V.^2);
CL=(L)./(dp.*c);
CD=(D)./(dp.*c);
disp(' ')
sprintf('Lift per unit span length L(N/m) is %f',L)
sprintf('Drag per unit span length D(N/m) is %f',D)
sprintf('Coefficient of lift is %f',CL)
sprintf('Coefficient of drag is %f',CD)
%% Save data to file
fileID = fopen(['coefficients_theoretical(Ma=',num2str(M),',alpha=',num2str(degrees),').txt'],'w');
fprintf(fileID,'Lift: %6.3f\n',CL);
fprintf(fileID,'Drag: %6.3f',CD);
fclose(fileID);
    
    
    


    

