classdef StefanSim2 < handle
    properties
        grid      % vector for the finite element grid 
        h %length of elements
        k %number of elements
        tmax %end time of simulation
        ph %physicalData
        T0 %dirichlet boundary function handle
        e %neumann boundary function handle
        Hini %inital enthalpy
        S %solution
        dt %timestep
        n %number of timesteps
        Mhd %dirichlet mass matrix
        R  %dirichlet vector
        lambda %for exact solution if known
        times %times where the solution is given
        Mh %mass matrix
        reg % use regularisation
        rfac % regularisation factor
        sys % number of linear systems solved
        fac %factor by which the standart timestep is increased/decreased
        diag %plot diagnostic variabes while running
        JPs % spares matrix template
        nr %first template vector for jacobi matrix
        fr %second template vector for jacobi matrix
        back % backward euler or crank nicolson?
    end
    methods
        function obj = StefanSim2(grid, tmax, ph, T0,e ,Hini,fac, diag,back ,varargin)

            % constructor for StefanSim2 class
            obj.back=back;
            obj.fac=fac;
            obj.grid = grid;
            obj.tmax = tmax;
            obj.diag=diag;


            if ~isa(ph, 'physicalData')
                error('ph must be a physicalData object');
            end
            obj.ph = ph;

            if ~isa(T0, 'function_handle')
                error('T0 must be a function handle');
            end
            obj.T0 = T0;

            if ~isa(e, 'function_handle')
                error('e must be a function handle');
            end
            obj.e = e;

            if ~isvector(Hini) || length(Hini) ~= length(grid)-1
                error('Hini must be a vector with length grid-1');
            end
            obj.Hini = Hini;
            obj.reg = false;
            if(length(varargin)>0)
            obj.rfac = varargin{1};
            end
        end

        function obj=initialize(obj)
            obj.k= length(obj.grid)-1;
            obj.h=obj.grid(2:end)-obj.grid(1:end-1);
            %obj.dt = obj.fac; % timestep
            obj.dt = min(obj.h)^2*obj.ph.c_fro(1,1)/obj.ph.k_fro(1,1)*obj.fac; % timestep

            obj.n = round(obj.tmax/obj.dt)+1;
            obj.S = zeros(obj.k, obj.n); % initialize solution
            obj.S(:,1)=obj.Hini;
            M = Mh(obj.h,true);
            obj.Mhd = M(2:end,2:end);
            obj.Mh = M(1:end,1:end);
            obj.R=M(2:end,1);
            obj.times=[(0:obj.dt:obj.tmax),obj.tmax];
            obj.JPs=speye(obj.k,obj.k);
            obj.nr=obj.ph.k_nor./obj.h;                                             
            obj.fr=obj.ph.k_fro./obj.h;  
            obj.findLam

        end

        function out = alpha(obj,x)
            % function that outputs the state of system, that is which
            % cells are frozen and which are unfrozen
            out=-(x<0)+(x>obj.ph.L);
        end

        function obj=simulate(obj)
            %function that performs the simulation

            n=obj.n; % numer of timesteps
            k= obj.k; 
            dt=obj.dt;
            H_c = obj.S(:,1); % current solution
            H_o=H_c;
            obj.sys=0; % total nuzmber of linear systems solved
            oldstate=2*ones(k,1);
            state=ones(k,1);
            x_o=ones(k,1);
            U=[];

            obj.JF(H_c)

            for i=1:n
                t_c = obj.dt*(i-1); % set current time

                if(obj.back)
                    Fn=@(x) FF(x,obj.T0(t_c+obj.dt),obj.ph,obj.h); % evaluate the function at the new time
                    C=obj.dt*(  ...  %Term C (Boundary condition)
                                       [zeros(k-1,1);obj.e(t_c+dt)] );  
    
                    Phi=@(x)obj.Mhd*(x-H_c)- ...                        %Term A (Mass Matrix)
                            obj.dt*(Fn(x))+ ...                %Term B (Stiffness part)
                            C;
                else
                    Fc=  FF(H_c,obj.T0(t_c),obj.ph,obj.h);         % evaluate the function at the old time
                    Fn=@(x) FF(x,obj.T0(t_c+obj.dt),obj.ph,obj.h); % evaluate the function at the new time
                    C=obj.dt/2*( [zeros(k-1,1);obj.e(t_c)] + ...  %Term C (Boundary condition)
                                 [zeros(k-1,1);obj.e(t_c+dt)] );  
                    
                    Phi=@(x)obj.Mhd*(x-H_c)- ...                        %Term A (Mass Matrix)
                            obj.dt/2*(Fn(x)+Fc)+ ...                %Term B (Stiffness part)
                            C;
                end

                res=1; % residual of newtons method
                x=H_c; % set initial guess to current timestep
               
     
                u=0; % number of newton iterations
                %b=Phi(x);   % rhight hand sight
                %K=zeros(k,1);
                %X=zeros(k,1);
                flag=false;
                while ~flag %res>10^-10*norm(x)    %%%%%   Newtons Method
                    u=u+1;       
                    % update the jacobi matrix
             
                    %obj.JF(x,k);                     
                    b=Phi(x);   % rhight hand sight

                    JP=updateJPhi(obj,k,state,oldstate);
                    %max(abs(eig(inv(obj.Mhd)*obj.JPs)) )
                    del=JP\b;   % solution of the linear system

                    alp=1;      % damping factor for netwons method
                    s=del;
                    
                    %Kratzneskow alg
                    
                    if(obj.ds(x-s,x)>0)
                        while(alp>10^-10) %%%%%% globalisation by sufficient decrease condition
                            
                            alp=alp/2;
                            if(obj.ds(x-s,x)>0 )
                                s=s-alp*del;
                     
                            else
                                s=s+alp*del;
                            end
                            norm(s,2);
                        end
                        while(obj.ds(x-s,x)==0&(alp>10^-14))
                            %alp=alp/2;
                            s=s+alp*del;
                        end
                        if(obj.ds(x-s,x)~=1)
                            s=s+(rand(size(x)) -ones(size(x))*0.5)*alp;
                        end

                    end
                    
                    %{
                    if(norm(obj.alpha(x-s)-obj.alpha(x),1)>1)
                        while(norm(obj.alpha(x-s)-obj.alpha(x),1)>1) %%%%%% globalisation by sufficient decrease condition
                            
                            alp=alp/2;
                            s=s-alp*del;
                        end
                        while(norm(obj.alpha(x-s)-obj.alpha(x),1)==0&(alp>10^-7))
                            alp=alp/2
                            s=s+alp*del;
                        end
                        if(norm(obj.alpha(x-s)-obj.alpha(x),1)>1)
                            s=s+rand(size(x))*10^-12;
                        end
                    end
                    %}
                    
                    
                    x=x-s;
                    %K=[K, obj.alpha(x)];
                    
                    %X=[X, b];
                
                    oldstate=state;
                    state=obj.alpha(x);
                    flag=obj.ds(x_o,x-s)==0&alp==1;

                    x_o=x;
                    
                    
                    if(obj.diag)
                    SaveAndPrintDiagnoseInfo(obj,u,x,res);
                    end
                end
                if(norm(Phi(x))>10^-6)
                    error("bigres")
                end
                %res=norm(Phi(x))
                %K=K(:,2:end-1);
                %{
                if( size(K,2)>8)
                    scatter(X(1,:),X(2,:), "x")
                    %error();
                end
                %}

                %X=X(:,2:end-1);

                if(obj.diag>0)
                    U=[U,u];
                    if(u>2)
                       u;
                       tu=obj.alpha(H_c)-obj.alpha(x);
                       tu(tu~=0);
                       
                    end
                end
                i
                obj.sys=obj.sys+u;
                H_o=H_c;
                H_c = x;

                obj.S(:,i+1)= x;
            end
            %U
            %mean(U)
            i;
            obj.sys;
        end

        function out=updateJPhi(obj,k,state,oldstate)

                    diff= find(state~=oldstate); %find columns of matrix that need to be uptaded

                    if(length(diff)>0)
                        
                        z=0;
                        q=0;
    
                        for(u=1:length(diff))
                            j=diff(u);
                            z=0;
                            q=0;
                            s=(state(j)==1)./obj.ph.c_nor(j)+ (state(j)==-1)./obj.ph.c_fro(j); %calculate scaling factor
                            
                            if(j<k)
                                z=(obj.nr(j+1).*(state(j)==1)+obj.fr(j+1).*(state(j)==-1)).*s; %calculate lower diag entry
                                obj.JPs(j+1,j)=z;
                            end
                            q=(obj.nr(j).*(state(j)==1)+obj.fr(j).*(state(j)==-1)).*s;%calculate upper diag entry
                            if(j>1)
                                
                                obj.JPs(j-1,j)=q;
                            end
                            obj.JPs(j,j)=-z-q;
                        end
                    end
                    if(obj.back)
                    out=obj.Mhd-obj.dt*obj.JPs;
                    else
                        out=obj.Mhd-0.5*obj.dt*obj.JPs;
                    end
                    %largEV=max(abs(eig(full(inv(obj.Mhd)*obj.JPs))))
                    %bigEl=max(max(abs(full(inv(obj.Mhd)*obj.JPs))))
                    %schur= sqrt(norm(full(inv(obj.Mhd)*obj.JPs),1)*norm(full(inv(obj.Mhd)*obj.JPs),"inf"))
                    %ershG= norm(full(inv(obj.Mhd)*obj.JPs),1)

                                            
        end

        function SaveAndPrintDiagnoseInfo(obj,u,x,res)

                    if(obj.diag>1)
                        u % print the current newton iteration
                        obj.alpha(x) % print the state of the system
                        res % print the residual
                    end

                    if(u>20)
                            %disp("newton iteration >20");
                    end
        end

        function JF(obj,x,k)

                al=obj.alpha(x);                                                    % system status
                s=(al==1)./obj.ph.c_nor+ (al==-1)./obj.ph.c_fro;      % column scaling factor

                nr=obj.ph.k_nor./obj.h;                                             % value for the vector r  for unfrozen cells
                fr=obj.ph.k_fro./obj.h;                                             % values for the vector r for frozen cells 

                ud=(nr.*(al==1)+fr.*(al==-1)).*s;                                      % upper diagonal 
                ld=(nr(2:end).*(al(1:end-1)==1)+fr(2:end).*(al(1:end-1)==-1)).*s(1:end-1);      % lower diagonal
                ld=[ld;0];
                d=-ud-ld;                                                       % diagonal

                obj.JPs=spdiags([ld,d,ud],[-1,0,1],obj.JPs);
                %obj.JPs(2:k,1:(k-1))=ld(1:end-1);% assemblation of out
                %obj.JPs(1:k-1,2:k)=ud(1:end);
                %obj.JPs(1:k,1:k)=d;

        end




        function flag = findLam(obj)
            obj.lambda = fzero(@(x)evallam2( x, obj.ph ,   obj.T0(0), TH(obj.Hini(1), obj.ph.c_fro(1,1),obj.ph.c_nor(1,1),obj.ph.L(1,1) ) )  ,1 );
        end


        function compareThistToExact(obj, x)
            dl=60*60*24;

            [~,node]=  min(abs(obj.grid-x)) ;
            if(obj.reg)
            nsol= THR( obj.S( node-1,:), obj.ph.c_fro(1,1), obj.ph.c_nor(1,1), obj.ph.L(1,1), obj.rfac);
            else
            nsol= TH( obj.S( node-1,:), obj.ph.c_fro(1,1), obj.ph.c_nor(1,1), obj.ph.L(1,1));

            end
            exsol= arrayfun(@(t) ...
                Tana(t,x,obj.ph, ...
                 obj.T0(0) ...
                ,TH(obj.Hini(1,1), obj.ph.c_fro(1,1),obj.ph.c_nor(1,1),obj.ph.L(1,1) ), ...
                obj.lambda), ...
                (0:length(nsol)-1)*obj.dt );
            Err=nsol-[exsol];
            %figure()
            if(obj.reg)

            plot((0:(length(Err)-1))*obj.dt/dl,Err,'color',"blue")
            else
            plot((0:(length(Err)-1))*obj.dt/dl,Err,'color',"red")
            end

            xlabel("Time in days")
            ylabel("T_{aprox} - T_{exact}")
            
        end


        function out=infnormError(obj, x)
            [~,node]=  min(abs(obj.grid-x)) ;
            if(obj.reg)
            nsol= THR( obj.S( node-1,:), obj.ph.c_fro(1,1), obj.ph.c_nor(1,1), obj.ph.L(1,1), obj.rfac);
            else
            nsol= TH( obj.S( node-1,:), obj.ph.c_fro(1,1), obj.ph.c_nor(1,1), obj.ph.L(1,1));

            end
            exsol= arrayfun(@(t) ...
                TanaE(t,x,obj.ph, ...
                 obj.T0(0) ...
                ,TH(obj.Hini(1,1), obj.ph.c_fro(1,1),obj.ph.c_nor(1,1),obj.ph.L(1,1) ), ...
                obj.lambda), ...
                (0:length(nsol)-1)*obj.dt );
            Err=nsol(2:end)-[exsol(2:end)];
            %figure()
           out=norm(Err,"inf")
            
        end
        
        function p=plotThist(obj, x,tscale,varargin)
            dl=60*60*24;


            [~,node]=  min(abs(obj.grid-x));

            if(obj.reg)
                nsol= THR( obj.S( node-1,:), obj.ph.c_fro(1,1), obj.ph.c_nor(1,1), obj.ph.L(1,1), obj.rfac) ;
            else
                nsol= TH( obj.S( node-1,:), obj.ph.c_fro(1,1), obj.ph.c_nor(1,1), obj.ph.L(1,1)) ;
            end
            

            if(isnumeric(tscale))
                x=obj.times/dl+tscale;
            else
                if(tscale=="days")
                    x=obj.times/dl;
                else
                    x=obj.times;
                end
            end

            z=min(length(x),length(nsol))
            p=plot(x(1:z),nsol(1:z), "LineWidth",1);

            if(isnumeric(tscale))
                
            else
            if(tscale=="days")
                xlabel("time in days")
            else
                xlabel("time in seconds")
            end
            end
            ylabel("T in °C")
            
        end

        function p=plotExact(obj,x,tscale)
        st=obj.times(1);
        en=obj.times(end);
        lintim=st:0.01:en;
            exsol= arrayfun(@(t) ...
                TanaE(t,x,obj.ph, ...
                 obj.T0(1) ...
                ,TH(obj.Hini(1,1), obj.ph.c_fro(1,1),obj.ph.c_nor(1,1),obj.ph.L(1,1) ), ...
                obj.lambda), ...
                lintim );            
            dl=60*60*24;

            if(tscale=="days")
                xlabel("time in days")
                p=plot(lintim/dl, exsol, LineWidth=0.8, Color="blue")

            else
                xlabel("time in seconds")
                p=plot(lintim, exsol, LineWidth=0.8)

            end
            ylabel("T in °c")
        end

        function out= ds(obj,a,b)
            out= norm(obj.alpha(a)-obj.alpha(b),1)*(norm(a-b,1)>10^-5);
        end

    end

end