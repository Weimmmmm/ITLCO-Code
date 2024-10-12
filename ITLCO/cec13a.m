function [ output ] = cec13a( x,D,func_num,O,M,shuff )
Pi=3.1415926535897932384626433832795029;
Exp=2.7182818284590452353602874713526625;
switch func_num
    case 1  %Sphere Function
        z=x-O;
        output=sum(z.^2);
    case 2  %Rotated High Conditioned Elliptic Function
        z=Tosz(M*(x-O),D);
        F2=0;
        for i=1:D
            f2=(10^6)^((i-1)/(D-1))*(z(i).^2);
            F2=F2+f2;
        end
        output=F2;
    case 3  %Rotated Bent Cigar Function
        z=M(D+1:2*D,:)*Tasy(M(1:D,:)*(x-O),0.5,D);
        output=z(1)^2+(10^6)*(sum(z.^2)-z(1)^2);
    case 4  %Rotated Discus Function
        z=Tosz(M*(x-O),D);
        output=(10^6)*z(1)^2+sum(z.^2)-z(1)^2;
    case 5  %Different Powers Function
        z=x-O;
        F5=0;
        for i=1:D
            f5=abs(z(i))^(2+4*(i-1)/(D-1));
            F5=F5+f5;
        end
        output=sqrt(F5);
    case 6  %Rotated Rosenbrock’s Function
        z=M*2.048*(x-O)/100+1;
        F6=0;
        for i=1:D-1
            f6=100*(z(i)^2-z(i+1))^2+(z(i)-1)^2;
            F6=F6+f6;
        end
        output=F6;
    case 7  %Rotated Schaffers F7 Function
        y=Dia(10,D)*M(D+1:2*D,:)*Tasy(M(1:D,:)*(x-O),0.5,D);
        for i=1:D-1
            z(i)=sqrt(y(i)^2+y(i+1)^2);
        end
        F7=0;
        for i=1:D-1
            f7=sqrt(z(i))+sqrt(z(i))*(sin(50*z(i)^0.2))^2;
            F7=F7+f7;
        end
        F7=(1/(D-1)*F7)^2;
        output=F7;
    case 8  %Rotated Ackley’s Function
        z=Dia(10,D)*M*Tasy(M*(x-O),0.5,D);
        ouput=-20*exp(-0.2*sqrt(sum(z.^2)/D))-exp(sum(cos(2*Pi.*z))/D)+20+Exp;
    case 9  %Rotated Weierstrass Function
        a=0.5;
        b=3;
        kmax=20;
        F91=zeros(1,D);
        F92=0;
        z=Dia(10,D)*M(D+1:2*D,:)*Tasy(M(1:D,:)*0.5*(x-O)/100,0.5,D);
        for i=1:D
            for k=0:kmax
                f91=a^k*cos(2*Pi*b^k*(z(i)+0.5));
                F91(i)=F91(i)+f91;
            end
        end
        for k=0:kmax
            f92=a^k*cos(2*Pi*b^k*0.5);
            F92=F92+f92;
        end
        output=sum(F91)-D*F92;
    case 10  %Rotated Griewank’s Function
        z=Dia(100,D)*M*600*(x-O)/100;
        F10=1;
        for i=1:D
            f10=cos(z(i)/sqrt(i));
            F10=F10*f10;
        end
        output=sum(z.^2)/4000-F10+1;
    case 11  %Rastrigin’s Function
        z=Dia(10,D)*Tasy(Tosz(5.12*(x-O)/100,D),0.2,D);
        F11=0;
        for i=1:D
            f11=z(i)^2-10*cos(2*Pi*z(i))+10;
            F11=F11+f11;
        end
        output=F11;
    case 12  %Rotated Rastrigin’s Function
        z=M(1:D,:)*Dia(10,D)*M(D+1:2*D,:)*Tasy(Tosz(M(1:D,:)*5.12*(x-O)/100,D),0.2,D);
        F12=0;
        for i=1:D
            f12=z(i)^2-10*cos(2*Pi*z(i))+10;
            F12=F12+f12;
        end
        output=F12;
    case 13  %Non-continuous Rotated Rastrigin’s Function
        xw=M*5.12*(x-O)/100;
        for i=1:D
            if abs(xw(i))<=0.5
                y(i)=xw(i);
            else
                y(i)=round(2*xw(i))/2;
            end
        end
        z=M*Dia(10,D)*M*Tasy(Tosz(y,D),0.2,D);
        F13=0;
        for i=1:D
            f13=z(i)^2-10*cos(2*Pi*z(i))+10;
            F13=F13+f13;
        end
        output=F13;
    case 14  %Schwefel’s Function
        g=zeros(1,D);
        z=Dia(10,D)*(1000*(x-O)/100)+420.9687462275036;
        for i=1:D
%              if abs(z(i))<=500%%%%%%%%%%按程序中编的
%                 g(i)=z(i)*sin(sqrt(abs(z(i))));
%             elseif z(i)>500
%                 g(i)=(500-mod(z(i),500))*sin(sqrt(500-mod(z(i),500)))-(z(i)-500)^2/(10000*D);%改了
%             else
%                 g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(500-mod(abs(z(i)),500)))-(z(i)+500)^2/(10000*D);%改了
%             end
%             if abs(z(i))<=500
%                 g(i)=z(i)*sin((abs(z(i)))^0.5);
%             elseif z(i)>500
%                 g(i)=(500-mod(z(i),500))*sin(sqrt(abs(500-mod(z(i),500))))-(z(i)-500)^2/(10000*D);
%             elseif z(i)<-500
%                 g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(abs(mod(abs(z(i)),500)-500)))-(z(i)+500)^2/(10000*D);
%             end
           if z(i)>500
                g(i)=(500-mod(z(i),500))*sin(sqrt(abs(500-mod(z(i),500))))-(z(i)-500)^2/(10000*D);%改了
           elseif z(i)<-500
                g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(abs(500-mod(abs(z(i)),500))))-(z(i)+500)^2/(10000*D);%改了
           else
                g(i)=z(i)*sin(sqrt(abs(z(i))));
           end
        end
        output=418.9828872724338*D-sum(g);
    case 15  %Rotated Schwefel’s Function
        g=zeros(1,D);
        z=Dia(10,D)*M*(1000*(x-O)/100)+420.9687462275036;
        for i=1:D
%              if abs(z(i))<=500%%%%%%%%%%%按程序编的
%                 g(i)=z(i)*sin(sqrt(abs(z(i))));
%             elseif z(i)>500
%                 g(i)=(500-mod(z(i),500))*sin(sqrt(500-mod(z(i),500)))-(z(i)-500)^2/(10000*D);%改了
%             else
%                 g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(500-mod(abs(z(i)),500)))-(z(i)+500)^2/(10000*D);%改了
%             end

%             if abs(z(i))<=500
%                 g(i)=z(i)*sin(abs(z(i))^0.5);
%             elseif z(i)>500
%                 g(i)=(500-mod(z(i),500))*sin(sqrt(abs(500-mod(z(i),500))))-(z(i)-500)^2/(10000*D);
%             elseif z(i)<-500
%                 g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(abs(mod(abs(z(i)),500)-500)))-(z(i)+500)^2/(10000*D);
%             end
           if z(i)>500
                g(i)=(500-mod(z(i),500))*sin(sqrt(500-mod(z(i),500)))-(z(i)-500)^2/(10000*D);%改了
           elseif z(i)<-500
                g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(500-mod(abs(z(i)),500)))-(z(i)+500)^2/(10000*D);%改了
           else
                g(i)=z(i)*sin(sqrt(abs(z(i))));
           end
        end
        output=418.9828872724338*D-sum(g);
    case 16  %Rotated Katsuura Function
        z=M*Dia(100,D)*(M*5*(x-O)/100);
        F161=1;
        for i=1:D
            F162=0;
            for j=1:32
                f162=abs(2^j*z(i)-round(2^j*z(i)))/2^j;
                F162=F162+f162;
            end
            f161=(1+i*F162)^(10/D^1.2);
            F161=F161*f161;
        end
        output=(10/D^2)*F161-10/D^2;
    case 17  %Lunacek bi-Rastrigin Function
        F171=0;
        F172=0;
        F173=0;
        miu0=2.5;
        d=1;
        s=1-(1/(2*sqrt(D+20)-8.2));
        miu1=-sqrt((miu0^2-d)/s);
        y=10*(x-O)/100;
        for i=1:D
            if O(i)<0
                xw(i)=-2*y(i)+miu0;
            elseif O(i)==0
                xw(i)=miu0;
            else
                xw(i)=2*y(i)+miu0;
            end
        end
        z=Dia(100,D)*(xw-miu0);
        for i=1:D
            f171=(xw(i)-miu0)^2;
            F171=F171+f171;
            f172=(xw(i)-miu1)^2;
            F172=F172+f172;
            f173=cos(2*Pi*z(i));
            F173=F173+f173;
        end
        output=min(F171,d*D+s*F172)+10*(D-F173);
    case 18  %Rotated Lunacek bi-Rastrigin Function
        F181=0;
        F182=0;
        F183=0;
        miu0=2.5;
        d=1;
        s=1-(1/(2*sqrt(D+20)-8.2));
        miu1=-sqrt((miu0^2-d)/s);
        y=10*(x-O)/100;
        for i=1:D
            if O(i)<0
                xw(i)=-2*y(i)+miu0;
            elseif O(i)==0
                xw(i)=miu0;
            else
                xw(i)=2*y(i)+miu0;
            end
        end
        z=M*Dia(100,D)*(M*(xw-miu0));
        for i=1:D
            f181=(xw(i)-miu0)^2;
            F181=F181+f181;
            f182=(xw(i)-miu1)^2;
            F182=F182+f182;
            f183=cos(2*Pi*z(i));
            F183=F183+f183;
        end
        output=min(F181,d*D+s*F182)+10*(D-F183);
    case 19  %Rotated Expanded Griewank’s plus Rosenbrock’s Function
        z=M*((x-O)*5/100)+1;
        BRF=zeros(1,D);
        for i=1:D-1
            BRF(i)=100*(z(i)^2-z(i+1))^2+(z(i)-1)^2;
        end
        BRF(D)=100*(z(D)^2-z(1))^2+(z(D)-1)^2;
        F191=0;
        F192=0;
        for i=1:D
            f191=(BRF(i)).^2;
            f192=cos(BRF(i));
            F191=F191+f191;
            F192=F192+f192;
        end
        output=F191/4000-F192+D; 
    case 20  %Rotated Expanded Scaffer’s F6 Function
        z=M(D+1:2*D,:)*Tasy(M(1:D,:)*(x-O),0.5,D);
        g=zeros(1,D);
        for i=1:D-1
            g(i)=0.5+((sin(sqrt(z(i)^2+z(i+1)^2)))^2-0.5)/(1+0.001*(z(i)^2+z(i+1)^2))^2;
        end
        g(D)=0.5+((sin(sqrt(z(D)^2+z(1)^2)))^2-0.5)/(1+0.001*(z(D)^2+z(1)^2))^2;
        output=sum(g);


end

