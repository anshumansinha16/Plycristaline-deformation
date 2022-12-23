 
function [] = mse658_cfd_code()
 



dt = 0.01;
 
n = 10; %input ('M');%scanf("Number of grid points in y", M); taken as 10.

p = 100; %input ('P'); % time step = 100 . i.e after 1 second.
p1 = 1000 ;%input ('P1'); % use gs to be 1000.

L= 1.63112803e-3; %0.001m
err = 0.001; 
dx = L/(n-1) ; % length wise
Xl=zeros(p,n);
 
 beta = 0.6;
b = 2.546e-10;
nu = 0.5;
alpha = 0.45;
mu= 0.025; %MPa
m= 0.39; % this variation is important , sensitivity of this term is maximum .thus we can iterate this between 0 to 0.50 , after this it goes down
n1= 1.37;
t = 100;
a = nu*mu*(b.^3);

Rs = 1e14;


 


R = zeros (1,n);
R0 = zeros(1,n);
V = zeros (1,n);
R1= zeros(1,n);
V(1,1) =0;
V(1,n) = 0;

 
flag1=0;
flag2=0;
flag3=0;

    
    for j = 1:n %(i=1;i<=M;i++)
        
        
        R0(1,j) = (10^5);
        R(1,j) = (10^5);
        R1(1,j)= (10^5);
        V(1,j) = (10^2);
    end
   
 V(1,1) =0;
V(1,n) = 0;   

iter2 =0;
 
while (iter2 <= p )%&& flag1==0)
    
    iter = 0;
    while( iter <= p1 )%&& flag3==0 )
        iter = iter +1;
        
        
 
            for j = 2:1:n-1
                
                x = L - (j-1)*dx ;
                %disp(x)
                V(1,j)=beta*b*nu*exp(((a)/(1.66e-24))*(abs(1-abs(((t*x)/(alpha*mu*b))).^m)).^n1);
                
                R(1,j) = R0(1,j) + (Rs*dt*0.45/(2*dx))*( -V(1,j+1) + V(1,j-1)) ;
                 
            end
               
    
        for j = 1:n %(i=1;i<=M;i++)
        if(abs(R(1,j)-R1(1,j))>=err)
          flag3 = flag3*0;
         else
          flag3 = flag3*1;
        end
        end
        
     R1 = R; 
     
    end
    
  
    %converging conditions
 
    for j = 1:n %(i=1;i<=M;i++)
        
        
      if(abs(R(1,j)-R0(1,j))>=err)
          flag1 = flag1*0;
      else
          flag1 = flag1*1;
      end
      
    end
    
    
    iter2 = iter2 +1;
    R0 = R;
   %disp(R)
    
   
       for j =1:n
       Xl(iter2,j) = R(1,j);
       end
   
   
end

if(flag3+flag2==0)
   disp('sucess') 
end


disp('we have chosen 10 grid points for discretisation  ');
disp('time-step taken as 0.01 second ');
disp('the following are the number of dislocations at each grid point for different time-steps ');


%f= Xl;
disp('number of dislocations')
disp(Xl*1e-6) %changing density into number



disp('we have chosen 10 grid points for discretisation  ');
disp('time-step taken as 0.01 second ');
disp('the above are the number of dislocations at each grid point for 100 time-steps ');


%filename='final.xlsx';
%xlswrite(filename,f)

end


