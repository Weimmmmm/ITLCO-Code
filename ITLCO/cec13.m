function [output] = cec13(x,D,func_num,O,M,shuff )
%cec2013 2019-7-8
%cec13命令 下列复制到自己的优化算法里（在第一次计算适应度值前，大循环外）
%对应变量要注意一致（如变量“维度”的名称为D，若为其他请修改成D）
%txt文件地址修改（如下列命令对应地址为F:\bullshit\bin\13input_data）
%种群矩阵为（D,N）,换句话说每列为一个个体（按照自己习惯自行修改）
% func_num=1;
% shift_filename='F:\bullshit\bin\13input_data\shift_data.txt';
% O=importdata(shift_filename);
% rotated_filename=['F:\bullshit\bin\13input_data\M_D',num2str(D),'.txt'];
% M=importdata(rotated_filename);
% shuff=0;
Pi=3.1415926535897932384626433832795029;
Exp=2.7182818284590452353602874713526625;
INF=1e99;
% EPS=1e-14;
o=O(:,1:D)';
switch func_num
    case 1  %Sphere Function
        z=x-o(:,1);
        output=sum(z.^2);
    case 2  %Rotated High Conditioned Elliptic Function
        z=Tosz(M(1:D,:)*(x-o(:,1)),D);
        F2=0;
        for i=1:D
            f2=(10^6)^((i-1)/(D-1))*(z(i)^2);
            F2=F2+f2;
        end
        output=F2;
    case 3  %Rotated Bent Cigar Function
        z=M(D+1:2*D,:)*Tasy(M(1:D,:)*(x-o(:,1)),0.5,D);
        output=z(1)^2+(10^6)*(sum(z.^2)-z(1)^2);
    case 4  %Rotated Discus Function
        z=Tosz(M(1:D,:)*(x-o(:,1)),D);
        output=(10^6)*z(1)^2+sum(z.^2)-z(1)^2;
    case 5  %Different Powers Function
        z=x-o(:,1);
        F5=0;
        for i=1:D
            f5=abs(z(i))^(2+4*(i-1)/(D-1));
            F5=F5+f5;
        end
        output=sqrt(F5);
    case 6  %Rotated Rosenbrock’s Function
        z=M(1:D,:)*2.048*(x-o(:,1))/100+1;
        F6=0;
        for i=1:D-1
            f6=100*(z(i)^2-z(i+1))^2+(z(i)-1)^2;
            F6=F6+f6;
        end
        output=F6;
    case 7  %Rotated Schaffers F7 Function
        y=Dia(10,D)*M(D+1:2*D,:)*Tasy(M(1:D,:)*(x-o(:,1)),0.5,D);
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
        z=Dia(10,D)*M(D+1:2*D,:)*Tasy(M(1:D,:)*(x-o(:,1)),0.5,D);
        output=-20*exp(-0.2*sqrt(sum(z.^2)/D))-exp(sum(cos(2*Pi.*z))/D)+20+Exp;
    case 9  %Rotated Weierstrass Function
        a=0.5;
        b=3;
        kmax=20;
        F91=zeros(1,D);
        F92=0;
        z=Dia(10,D)*M(D+1:2*D,:)*Tasy(M(1:D,:)*0.5*(x-o(:,1))/100,0.5,D);
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
        z=Dia(100,D)*M(1:D,:)*600*(x-o(:,1))/100;
        F10=1;
        for i=1:D
            f10=cos(z(i)/sqrt(i));
            F10=F10*f10;
        end
        output=sum(z.^2)/4000-F10+1;
    case 11  %Rastrigin’s Function
        z=Dia(10,D)*Tasy(Tosz(5.12*(x-o(:,1))/100,D),0.2,D);
        F11=0;
        for i=1:D
            f11=z(i)^2-10*cos(2*Pi*z(i))+10;
            F11=F11+f11;
        end
        output=F11;
    case 12  %Rotated Rastrigin’s Function
        z=M(1:D,:)*Dia(10,D)*M(D+1:2*D,:)*Tasy(Tosz(M(1:D,:)*5.12*(x-o(:,1))/100,D),0.2,D);
        F12=0;
        for i=1:D
            f12=z(i)^2-10*cos(2*Pi*z(i))+10;
            F12=F12+f12;
        end
        output=F12;
    case 13  %Non-continuous Rotated Rastrigin’s Function
        xw=M(1:D,:)*5.12*(x-o(:,1))/100;
        for i=1:D
            if abs(xw(i))<=0.5
                y(i)=xw(i);
            else
                y(i)=round(2*xw(i))/2;
            end
        end
        z=M(1:D,:)*Dia(10,D)*M(D+1:2*D,:)*Tasy(Tosz(y,D),0.2,D);
        F13=0;
        for i=1:D
            f13=z(i)^2-10*cos(2*Pi*z(i))+10;
            F13=F13+f13;
        end
        output=F13;
    case 14  %Schwefel’s Function
        g=zeros(1,D);
        z=Dia(10,D)*(1000*(x-o(:,1))/100)+4.209687462275036e+002;
        for i=1:D
%              if abs(z(i))<=500%%%%%%%%%%%按程序中编的
%                 g(i)=z(i)*sin(sqrt(abs(z(i))));
%             elseif z(i)>500
%                 g(i)=(500-mod(z(i),500))*sin(sqrt(500-mod(z(i),500)))-(z(i)-500)^2/(10000*D);%改了
%             else
%                 g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(500-mod(abs(z(i)),500)))-(z(i)+500)^2/(10000*D);%改了
%             end

%             if abs(z(i))<=500%%%%%%%%%%%%%按PDF中编的
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
        output=4.189828872724338e2*D-sum(g);
   case 15  %Rotated Schwefel’s Function
        g=zeros(D,1);
        z=Dia(10,D)*M(1:D,:)*(1000*(x-o(:,1))/100)+4.209687462275036e2;
        for i=1:D
%              if abs(z(i))<=500%%%%%%%%%按PDF中编的
%                 g(i)=z(i)*sin((abs(z(i)))^0.5);
%             elseif z(i)>500
%                 g(i)=(500-mod(z(i),500))*sin(sqrt(abs(500-mod(z(i),500))))-(z(i)-500)^2/(10000*D);
%             elseif z(i)<-500
%                 g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(abs(mod(abs(z(i)),500)-500)))-(z(i)+500)^2/(10000*D);
%             end

%             if abs(z(i))<=500%%%%%%%%%%按程序中编的
%                 g(i)=z(i)*sin(sqrt(abs(z(i))));
%             elseif z(i)>500
%                 g(i)=(500-mod(z(i),500))*sin(sqrt(500-mod(z(i),500)))-(z(i)-500)^2/(10000*D);%改了
%             else
%                 g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(500-mod(abs(z(i)),500)))-(z(i)+500)^2/(10000*D);%改了
%             end
           
           if z(i)>500
                g(i)=(500-mod(z(i),500))*sin(sqrt(abs(500-mod(z(i),500))))-(z(i)-500)^2/(10000*D);%改了
           elseif z(i)<-500
                g(i)=(mod(abs(z(i)),500)-500)*sin(sqrt(abs(500-mod(abs(z(i)),500))))-(z(i)+500)^2/(10000*D);%改了
           else
                g(i)=z(i)*sin(sqrt(abs(z(i))));
           end
        end
        output=4.189828872724338e2*D-sum(g);
    case 16  %Rotated Katsuura Function
        z=M(D+1:2*D,:)*Dia(100,D)*(M(1:D,:)*5*(x-o(:,1))/100);
        F161=1;
        for i=1:D
            F162=0;
            for j=1:32
                f162=abs((2^j)*z(i)-round((2^j)*z(i)))/2^j;
                F162=F162+f162;
            end
            f161=(1+i*F162)^(10/(D^1.2));
            F161=F161*f161;
        end
        output=(10/(D^2))*F161-10/(D^2);
    case 17  %Lunacek bi-Rastrigin Function
        F171=0;
        F172=0;
        F173=0;
        miu0=2.5;
        d=1;
        s=1-(1/(2*sqrt(D+20)-8.2));
        miu1=-sqrt((miu0^2-d)/s);
        y=10*(x-o(:,1))/100;
        xw=zeros(D,1);
        for i=1:D
            if o(i,7)<0
                xw(i,1)=-2*y(i)+miu0;
            else
                xw(i,1)=2*y(i)+miu0;
            end
        end
        z=Dia(100,D)*(xw-miu0);
        for i=1:D
            f171=(xw(i,1)-miu0)^2;
            F171=F171+f171;
            f172=(xw(i,1)-miu1)^2;
            F172=F172+f172;
            f173=cos(2*Pi*z(i));
            F173=F173+f173;
        end
%         F172=d*D+s*F172;
%         if F171>F172
%            F170=F172;
%         else F170=F171;
%         end
        output=min(F171,d*D+s*F172)+10*(D-F173);
    case 18  %Rotated Lunacek bi-Rastrigin Function
        F181=0;
        F182=0;
        F183=0;
        miu0=2.5;
        d=1;
        s=1-(1/(2*sqrt(D+20)-8.2));
        miu1=-sqrt((miu0^2-d)/s);
        y=10*(x-o(:,1))/100;
        xw=zeros(D,1);
        for i=1:D
            if o(i,8)<0
                xw(i,1)=-2*y(i)+miu0;
            elseif o(i,8)==0
                xw(i,1)=miu0;
            else
                xw(i,1)=2*y(i)+miu0;
            end
        end
        z=M(D+1:2*D,:)*(Dia(100,D)*(M(1:D,:)*(xw-miu0)));
        for i=1:D
            f181=(xw(i,1)-miu0)^2;
            F181=F181+f181;
            f182=(xw(i,1)-miu1)^2;
            F182=F182+f182;
            f183=cos(2*Pi*z(i));
            F183=F183+f183;
        end
%         F182=d*D+s*F182;
%         if F181>F182
%             F180=F182;
%         else F180=F181;
%         end
        output=min(F181,d*D+s*F182)+10*(D-F183);
    case 19  %Rotated Expanded Griewank’s plus Rosenbrock’s Function
        z=M(1:D,:)*((x-o(:,1))*5/100)+1;
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
       z=M(D+1:2*D,:)*Tasy(M(1:D,:)*(x-o(:,1)),0.5,D);
        g=zeros(1,D);
        for i=1:D-1
            g(i)=0.5+((sin(sqrt(z(i)^2+z(i+1)^2)))^2-0.5)/(1+0.001*(z(i)^2+z(i+1)^2))^2;
        end
        g(D)=0.5+((sin(sqrt(z(D)^2+z(1)^2)))^2-0.5)/(1+0.001*(z(D)^2+z(1)^2))^2;
        output=sum(g);
    case 21  %Composition Function 1
        N=5;
        output=0;
        derta=[10 20 30 40 50];
        lamuda=[1 1e-6 1e-26 1e-6 0.1];
        bias=[0 100 200 300 400];
         w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,6,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,5,O(2,1:D)',M(D+1:2*D,:),shuff);
        g(3)=cec13a(x,D,3,O(3,1:D)',M(2*D+1:4*D,:),shuff);
        g(4)=cec13a(x,D,4,O(4,1:D)',M(4*D+1:5*D,:),shuff);
        g(5)=cec13a(x,D,1,O(5,1:D)',M(5*D+1:6*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;
    case 22  %Composition Function 2
         N=3;
        output=0;
        derta=[20 20 20];
        lamuda=[1 1 1];
        bias=[0 100 200];
         w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,14,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,14,O(2,1:D)',M(D+1:2*D,:),shuff);
        g(3)=cec13a(x,D,14,O(3,1:D)',M(2*D+1:3*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;
    case 23  %Composition Function 3
        N=3;
        output=0;
        derta=[20 20 20];
        lamuda=[1 1 1];
        bias=[0 100 200];
        w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,15,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,15,O(2,1:D)',M(D+1:2*D,:),shuff);
        g(3)=cec13a(x,D,15,O(3,1:D)',M(2*D+1:3*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;
    case 24  %Composition Function 4
        N=3;
        output=0;
        derta=[20 20 20];
        lamuda=[0.25 1 2.5];
        bias=[0 100 200];
        w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,15,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,12,O(2,1:D)',M(D+1:3*D,:),shuff);
        g(3)=cec13a(x,D,9,O(3,1:D)',M(2*D+1:4*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;
    case 25  %Composition Function 5
        N=3;
        output=0;
        derta=[10 30 50];
        lamuda=[0.25 1 2.5];
        bias=[0 100 200];
         w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,15,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,12,O(2,1:D)',M(D+1:3*D,:),shuff);
        g(3)=cec13a(x,D,9,O(3,1:D)',M(2*D+1:4*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;
    case 26  %Composition Function 6
         N=5;
        output=0;
        derta=[10 10 10 10 10];
        lamuda=[0.25 1 1e-7 2.5 10];
        bias=[0 100 200 300 400];
        w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,15,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,12,O(2,1:D)',M(D+1:3*D,:),shuff);
        g(3)=cec13a(x,D,2,O(3,1:D)',M(2*D+1:3*D,:),shuff);
        g(4)=cec13a(x,D,9,O(4,1:D)',M(3*D+1:5*D,:),shuff);
        g(5)=cec13a(x,D,10,O(5,1:D)',M(5*D+1:6*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;
    case 27  %Composition Function 7
        N=5;
        output=0;
        derta=[10 10 10 20 20];
        lamuda=[100 10 2.5 25 0.1];
        bias=[0 100 200 300 400];
         w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,10,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,12,O(2,1:D)',M(D+1:3*D,:),shuff);
        g(3)=cec13a(x,D,15,O(3,1:D)',M(2*D+1:3*D,:),shuff);
        g(4)=cec13a(x,D,9,O(4,1:D)',M(3*D+1:5*D,:),shuff);
        g(5)=cec13a(x,D,1,O(5,1:D)',M(5*D+1:6*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;
    case 28  %Composition Function 8
        N=5;
        output=0;
        derta=[10 20 30 40 50];
        lamuda=[2.5 2.5e-3 2.5 5e-4 0.1];
        bias=[0 100 200 300 400];
        w_max=0;
        for i=1:N
            Ww(i)=sum((x-o(:,i)).^2);
            if Ww(i)==0
            W(i)=INF;
            else
            W(i)=(1/sqrt(Ww(i)))*exp(-Ww(i)/(2*D*derta(i)^2));
            end
            if W(i)>w_max
                w_max=W(i);
            end
        end
        for i=1:N
            SUM=sum(W);
        end
        if w_max==0
            for i=1:N
                W(i)=1;
                SUM=N;
            end
        end
        for i=1:N
            w(i)=W(i)/SUM;
        end
        g(1)=cec13a(x,D,19,O(1,1:D)',M(1:D,:),shuff);
        g(2)=cec13a(x,D,7,O(2,1:D)',M(D+1:3*D,:),shuff);
        g(3)=cec13a(x,D,15,O(3,1:D)',M(2*D+1:3*D,:),shuff);
        g(4)=cec13a(x,D,20,O(4,1:D)',M(3*D+1:5*D,:),shuff);
        g(5)=cec13a(x,D,1,O(5,1:D)',M(5*D+1:6*D,:),shuff);
        for i=1:N
            f=w(i)*(lamuda(i)*g(i)+bias(i));
            output=output+f;
        end
        output=output;

end

