function [t,y_out] = ND_integrator(fn,t,y0,K,M_inv,omega,integrator)%,options)

%% Helper functions
    function LMM(y0,t)
        h = t(2)-t(1);
        if(sum(abs(t(2:2:end)-t(1:2:end-1)-h)>1e-6)>0)
            error('For LMM, time vector must have equidistant nodes!');
        end
        err=1e6;
        iter_total=0;
        while(err>err_tol)
            iter_total=iter_total+1;
            t_init_temp = t(1)+(0:3).*h;
            y_in=y_out;
            rkf45(y0,t_init_temp);
            y_init = y_out(1:4,:)';
            y_out = y_in;
            y_out(:,(t_init_temp-t(1:4)<1e-6))=y_init(:,(t_init_temp-t(1:4)<1e-6));
            step_total = round((t(end)-t(1))/h);
            t_curr = t_init_temp;
            y_curr = y_init;
            while(t_curr+h < t(end))
                % initialize space and time for next iteration
                y_old = y_curr;
                y_curr(:,1:3) = y_old(:,2:4);
                t_old = t_curr;
                t_curr=t_curr+h;
                % Adams-Bashforth Predictor
                p=y_old(:,4)+h/24*(55*fn(t_old(4),y_old(:,4),K,M_inv,omega)-...
                    59*fn(t_old(3),y_old(:,3),K,M_inv,omega)+...
                    37*fn(t_old(2),y_old(:,2),K,M_inv,omega)-...
                    9*fn(t_old(1),y_old(:,1),K,M_inv,omega));
                % Adams-Moulton Corrector
                y_curr(:,4)=y_old(:,4)+h/24*(9*fn(t_curr(4),p,K,M_inv,omega)+...
                    19*fn(t_old(4),y_old(:,4),K,M_inv,omega)-...
                    5*fn(t_old(3),y_old(:,3),K,M_inv,omega)+...
                    fn(t_old(2),y_old(:,2),K,M_inv,omega));
                err = (19/720)*(sum((y_curr(:,4)-p).^2)/(norm(y_curr(:,4))+sqrt(eps)));
                if(err>err_tol)
                    break;
                end
                % print stuff - for troubleshooting
                if(mod(iter_total,1e2)==0)
                    fprintf(1,'h: %.5e, R: %.5e\n',dh,R);
                    fprintf(1,'step: %d, time: %.2e s\n',rk_idx,t_temp);
                    fprintf(1,'while_loop iters: %d\n',iters);
                end
                if(sum(t==t_curr(4))>0)
                    index = find(t==t_curr(4),1,'first');
                    y_out(:,index) = y_curr(:,4);
                end
            end
            h=h/2;
            if(h<step_tol)
                fprintf(1,'h: %.5e, Error: %.5e\n',h,err);
                fprintf(1,'while_loop iters: %d\n',iter_total);
                error('h is too small!');
            end
        end
        fprintf(1,'Total LMM while loops: %d\n',iter_total);
        fprintf(1,'Total LMM steps: %d\n',step_total);
        y_out = y_out';
    end

    function g = LMMMinFun(x0,t,y0,h)
%         g = x0 - ((1/25).*(48.*y0(4,:)-36.*y0(3,:)...
%                 +16.*y0(2,:)-3.*y0(1,:)+...
%                 12.*h.*fn(t,x0',K,M_inv,omega)'));
%         g=g*g';
        g=x0-y0(:,4)+h/24*(9*fn(t(5),p,K,M_inv,omega)+...
            19*fn(t(4),y0(:,4),K,M_inv,omega)-...
            5*fn(t(3),y0(:,3),K,M_inv,omega)+...
            fn(t(2),y0(:,2),K,M_inv,omega));
    end

    function bd4(y0,t)
        rkf45(y0,t(1:4));
        bd_t = t;
        for bd_idx=5:length(bd_t)
            h = bd_t(bd_idx)-bd_t(bd_idx-1);
%             y_out(:,bd_idx) = (1/25).*(48.*y(:,bd_idx-1)-36.*y(:,bd_idx-2)...
%                 +16.*y(:,bd_idx-3)-3.*y(:,bd_idx-4)+...
%                 12.*h.*fn(bd_t,y_out(:,bd_idx),K,M_inv,omega));
            y_out(bd_idx,:) = fminsearch(@(x0) bdMinFun(x0,t,...
                y_out(bd_idx-4:bd_idx-1,:),h),y_out(bd_idx-1,:))';
%             fprintf(1,'Next! %d\n',bd_idx);
        end
    end
    
    function g = bdMinFun(x0,t,y0,h)
        g = x0 - ((1/25).*(48.*y0(4,:)-36.*y0(3,:)...
                +16.*y0(2,:)-3.*y0(1,:)+...
                12.*h.*fn(t,x0',K,M_inv,omega)'));
        g=g*g';
    end
  
    function rk4(y0,t)
        % fixed step rk4
        % input requires time vector, whose intervals may be non-fixed
        rk_t=t;
        y_out(1,:)=y0;
        for rk_idx=1:(length(rk_t)-1)
            % initialize step size
            h=t(rk_idx+1)-t(rk_idx);
            dh = h;
            % initialize temp vars (less typing)
            y_temp = y_out(rk_idx,:);
            t_temp = rk_t(rk_idx+1);
            % calculate 5 stage rk for fourth order accuracy
            k1=dh.*fn(t_temp,y_temp',K,M_inv,omega)';
            k2=dh.*fn(t_temp+dh/4,(y_temp+k1./4)',K,M_inv,omega)';
            k3=dh.*fn(t_temp+3*dh/8,(y_temp+(3/32).*k1+(9/32).*k2)'...
                ,K,M_inv,omega)';
            k4=dh.*fn(t_temp+12*dh/13,(y_temp+(1932/2197).*k1-...
                (7200/2197).*k2+(7296/2197).*k3)',K,M_inv,omega)';
            k5=dh.*fn(t_temp+dh,(y_temp+(439/216).*k1-(8).*k2+...
                (3680/513).*k3-(845/4104).*k4)',K,M_inv,omega)';
            % calculate new step
            y_new = y_temp+(25/216).*k1+(1408/2565).*k3+...
                (2194/4104).*k4-(1/5).*k5;
            % set next value
            y_out(rk_idx+1,:) = y_new;
        end
    end

    function rkf45(y0,t)
        iter_total=0;
        step_total=0;
        rk_t=t;
        y_out(:,1)=y0;
        h=t(2)-t(1);
        dh = h*0.9;% try almost full step
        for rk_idx=1:(length(rk_t)-1)
            % initialize space and time for next iteration
            y_temp = y_out(:,rk_idx);
            t_temp = rk_t(rk_idx);
            % count while loop iterations
            iters=0;
            % count sub-steps taken
            steps=0;
            % while we haven't gone to the next time step...
            while(1==1)
                % count loop iterations
                iters = iters+1;
                % setup RK stages
                k1=dh.*fn(t_temp,y_temp,K,M_inv,omega);
                k2=dh.*fn(t_temp+dh/4,(y_temp+k1./4),K,M_inv,omega);
                k3=dh.*fn(t_temp+3*dh/8,(y_temp+(3/32).*k1+(9/32).*k2)...
                    ,K,M_inv,omega);
                k4=dh.*fn(t_temp+12*dh/13,(y_temp+(1932/2197).*k1-...
                    (7200/2197).*k2+(7296/2197).*k3),K,M_inv,omega);
                k5=dh.*fn(t_temp+dh,(y_temp+(439/216).*k1-(8).*k2+...
                    (3680/513).*k3-(845/4104).*k4),K,M_inv,omega);
                k6=dh.*fn(t_temp+dh/2,(y_temp-(8/27).*k1+(2).*k2-...
                    (3544/2565).*k3+(1859/4104).*k4-(11/40).*k5),K,M_inv,omega);
                % calculate RK steps
                y_new = y_temp+(25/216).*k1+(1408/2565).*k3+...
                    (2194/4104).*k4-(1/5).*k5;
                y_hat = y_temp+(16/135).*k1+(6656/12825).*k3+...
                    (28561/56430).*k4-(9/50).*k5+(2/55).*k6;
                % calculate error (LTE)
                R = (1/dh)*sum((y_hat-y_new).^2);
                % rescaling factor
                delta = 0.84*(err_tol/R)^0.25;
                if(R<=err_tol)% if error is under control
                    steps=steps+1;
                    if(t_temp==rk_t(rk_idx+1))
                        y_temp = y_new;
                        break;
                    elseif(t_temp+dh>rk_t(rk_idx+1))
                        dh = rk_t(rk_idx+1)-t_temp;
                    end
                    y_temp = y_new;
                    t_temp = t_temp+dh;
                % say something if we are failing to converge
                elseif(dh<step_tol)
                    fprintf(1,'h: %.5e, R: %.5e\n',dh,R);
                    fprintf(1,'step: %d, time: %.2e s\n',rk_idx,t_temp);
                    fprintf(1,'while_loop iters: %d\n',iters);
                    error('h is too small!');
                end
                % print stuff - for troubleshooting
                if(mod(iters,1e4)==0)
                    fprintf(1,'h: %.5e, R: %.5e\n',dh,R);
                    fprintf(1,'step: %d, time: %.2e s\n',rk_idx,t_temp);
                    fprintf(1,'while_loop iters: %d\n',iters);
                end
                % rescale time step
                dh = delta*dh;
            end
            y_out(:,rk_idx+1) = y_temp;
            iter_total=iter_total+iters;
            step_total=step_total+steps;
        end
        fprintf(1,'Total rkf45 while loops: %d\n',iter_total);
        fprintf(1,'Total rkf45 steps: %d\n',step_total);
        y_out = y_out';
    end

%% Setup
fn = str2func(fn);
n = length(t);
if(n==2)
    n=1000;
    t = linspace(t(1),t(2),n);
elseif(n<2)
    error('Time vector should have at least two points!');
end
m = length(y0);
y_out = zeros(m,n);

step_tol = 1e-8;
err_tol = 1e-6;

integrator = str2func(['ND_integrator/',integrator]);

%% Main call
% use Runge-Kutta method; rkf45 is the adaptive routine, while rk4 is fixed 
integrator(y0,t);

end