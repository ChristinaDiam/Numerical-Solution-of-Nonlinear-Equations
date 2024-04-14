%% Diamanti Xristina 1115201800046

%% 1.1

f1=inline('(x+1).^3.*(x-2)');           %f1
df1=inline('(x+1).^2.*(4*x-5)');    %paragogos f1
a1_1=0;
b1_1=5;

NRresult1=syndyasmos_D_NR(f1,df1,a1_1,b1_1) %thetiki riza


%% 1.1

f1=inline('(x+1).^3.*(x-2)');           %f1
df1=inline('(x+1).^2.*(4*x-5)');    %paragogos f1
a1_2=-3;
b1_2=0;

NRresult2=syndyasmos_D_NR(f1,df1,a1_2,b1_2) %arnitiki riza


%% 1.1

f2=inline('exp(x)-x.^2-2');           %f2
df2=inline('exp(x)-2*x');    %paragogos f2
a2=0;
b2=5;

r=syndyasmos_D_NR(f2,df2,a2,b2)


%% 1.3 

[en1, posotita1_1, posotita2_1] = siglisi(NRresult1,2) %gia f1 (arnitiki riza)


%% 1.3

[en2 posotita1_2 posotita2_2] = siglisi(NRresult2,2) %gia f1 (thetiki riza)


%% 1.3

[en3 posotita1_3 posotita2_3] = siglisi(r,2) %gia f2


%% 1.5

f1=inline('(x+1).^3.*(x-2)');           %f1
a1_1=0;
b1_1=5;

DTresult1=syndyasmos_D_T(f1,a1_1,b1_1) %thetiki riza


%% 1.5

f1=inline('(x+1).^3.*(x-2)');           %f1
a1_2=-3;
b1_2=0;

DTresult2=syndyasmos_D_T(f1,a1_2,b1_2)  %arnitiki riza


%% 1.5

f2=inline('exp(x)-x.^2-2');           %f2
a2=0;
b2=5;

rDT=syndyasmos_D_T(f2,a2,b2) 


%% 1.5

[en1_5 posotita1_1_5 posotita2_1_5] = siglisi(DTresult1,2) %gia f1 (thetiki riza)


%% 1.5

[en2_5 posotita1_2_5 posotita2_2_5] = siglisi(DTresult2,1.6) %gia f1 (arnitiki riza)


%% 1.5

[en3_5 posotita1_3_5 posotita2_3_5] = siglisi(rDT,1.6) %gia f2


%% functions

function out_nr=syndyasmos_D_NR(f,df,a,b)
    ed=0.5*10^(-2);
    eNR=0.5*10^(-6);
   
    out_d=bisect_m(f, a, b, ed, 50)    %dixotomisi
    
    %telefteo diastima ap tin dixotomisi
    a_end=out_d(end,2);
    b_end=out_d(end,3);
    
    m=(b_end+a_end)/2;      %meso diastimatos
   
    out_nr=rf_newton2(f, df, m, eNR, 50) %NR

end


function [en, posotita1, posotita2] = siglisi(NRresult,p)
    riza=NRresult(end,2);
    
    for n=1:height(NRresult)
        en(n)=abs(NRresult(n,2)-riza); %apolito sfalma (gnosti riza)
        
        if n<height(NRresult)
            posotita1(n)=abs(NRresult(n+1,2)-riza)/abs(en(n))^p;
        end
        
        if (n>=2 && n<height(NRresult))
            posotita2(n)=abs(NRresult(n+1,2)-NRresult(n,2))/abs(NRresult(n,2)-NRresult(n-1,2))^p;
        end
    end
end


function out_temn=syndyasmos_D_T(f,a,b)
    ed=0.5*10^(-2);
    eT=0.5*10^(-6);
   
    out_d=bisect_m(f, a, b, ed, 50)    %dixotomisi
    
    %telefteo diastima ap tin dixotomisi
    a_end=out_d(end,2);
    b_end=out_d(end,3);
   
    out_temn=temnousa(f, b_end, a_end, eT, 50)

end